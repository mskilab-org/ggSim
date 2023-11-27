#' @import skitools
#' @import skidb
#' @import data.table
#' @import gGnome
#' @import gUtils
#' @import bamUtils

#' @name .simcov
#' @title 
#'
#' @description
#'
#' @param gr 
#' @param binsize 
#' @param bins 
#' @param bias 
#' @param diploid 
#' @param readsize 
#' @param basecov 
#' @param overdispersion 
#' @param purity 
#' @param poisson 
#' @param normalize 
.simcov = function(gr, ## needs field 'cn'
                   binsize = 1e3,
                   bins = gr.tile(seqlengths(gr), binsize),
                   bias = NULL, ## granges of biases for given regions, first field is assumed to be the bias amount
                   diploid = TRUE,
                   readsize = 150,
                   basecov = 60,
                   overdispersion = NA,
                   purity = 0.95,
                   poisson = TRUE,
                   normalize = TRUE)
{
  if (inherits(gr, 'gGraph'))
    gr = gr$nodes$gr
  
  binsize = unique(width(bins))[1]
  
  if (diploid)
    ncn = 2
  else
    ncn = 1
  
  bincov = (basecov*(readsize+binsize-1))/readsize ## estimate num reads per bin
  
  bins$cn = (bins %$% gr[, 'cn'])$cn
  bins$cn.rel = (purity*bins$cn + ncn*(1-purity))/(purity*mean(bins$cn, na.rm = TRUE)+ ncn*(1-purity))
  
  if (poisson) {
    if (is.na(overdispersion))
      bins$cov = rpois(length(bins), bins$cn.rel*bincov)
    else
      bins$cov = rnbinom(length(bins), mu = bins$cn.rel*bincov,
                         size = 1/overdispersion)
  } else {
    bins$cov = bins$cn.rel * bincov
  }
  
  if (normalize)
    bins$cov = bins$cov/mean(bins$cov, na.rm = TRUE)
  
  
  if (!is.null(bias))
  {
    values(bias)[[1]] = values(bias)[[1]]/mean(values(bias)[[1]], na.rm = TRUE)
    bins$cov.og = bins$cov
    ##bins$cov = bins$cov * values(bins[, c()] %$% bias[, 1])[[1]] ## made this immune to NA for targeted seq (Addy)
    bins$cov = bins$cov * values(gr.val(bins[, c()], bias[, 1], na.rm = T, val = names(values(bias[, 1]))))[[1]]
  }
  
  return(bins)
}

#' @name .separate
#'
#' @param dt 
#' @param old 
#' @param new 
#' @param sep 
#' @param perl 
#'
#' @return
#' @export
#'
#' @examples
.separate = function(dt, old, new, sep, perl = FALSE)
{
  newdt = dt[[old]] %>% as.character %>% strsplit(split = sep, perl = perl) %>% do.call(rbind, .) %>% as.data.table %>% setnames(new)
  dt[[old]] = NULL
  dt = cbind(dt, newdt)
  return(dt)
}

#' @name ggsim
#'
#' @param junctions Path to junction file e.g. bedpe, vcf, rds
#' @param vcf Phased VCF of germline hets
#' @param bias rds of binned read depth bias for tumor sample e.g. read depth for a random normal sample
#' @param nbias rds of binned read depth bias for normal sample e.g. read depth for a random normal sample
#' @param snps Optional comprehensive VCF of reference snps e.g. hapmap
#' @param unmappable Optional .rds of GRanges of CN unmappable regions
#' @param coverage target tumor base coverage
#' @param ncoverage target normal base coverage
#' @param alpha Target purity
#' @param poisson add poisson noise
#' @param tau Target ploidy
#' @param numbreaks number of additional breaks to provide in cn unmappable regions
#' @param width Binwidth for binned read depth
#' @param cnloh add a cnloh edge
#' @param outdir output directory
#' @param libdir Directory containing this R file
#'
#' @return
#' @export
#'
#' @examples
ggsim = function(junctions,
                 vcf,
                 bias,
                 nbias,
                 snps,
                 unmappable,
                 coverage,
                 ncoverage,
                 alpha,
                 poisson,
                 tau,
                 numbreaks,
                 width,
                 cnloh,
                 outdir,
                 libdir)
{
  if (is.null(libdir) | (is.null(junctions)) | is.null(vcf) | is.null(bias)| is.null(nbias))
  stop()
  
  setDTthreads(1)
  
  system(paste('mkdir -p',  outdir))
  
  if (!is.null(unmappable))
    CNun = unmappable %>% readRDS
  
  message('Loading phased SNPs')
  snps  = skidb::read_vcf(vcf, geno = TRUE)
  fn = rev(names(values(snps)))[1]
  snps = snps[lengths(snps$ALT)==1] ## remove multiallelic
  snps$ALT = unstrsplit(snps$ALT)
  snps = snps[nchar(snps$REF)==1 & nchar(snps$ALT) == 1]
  snps = snps[, c('REF', 'ALT', fn)] %>% as.data.frame %>% .separate(fn, c('A', 'B'), '\\|') %>% dt2gr %>% gr.sub
  if (grepl("rds$", junctions)) {
    junctions = jJ(readRDS(junctions))
  } else {
    junctions = jJ(junctions)
  }
  snps$het = snps$A != snps$B
  snpmat = values(snps)[, c('REF', 'ALT')] %>% as.matrix
  snps$A = snpmat[cbind(1:length(snps), as.numeric(snps$A)+1)] ## convert integer phases to bases
  snps$B = snpmat[cbind(1:length(snps), as.numeric(snps$B)+1)]
  message('ingested junctions and phased SNPs')
  
  ## additional snps that are homozogyous REF in the reference sample may not be
  ## present in the phased VCF, so we can pull these from a reference DB e.g. hapmap
  if (!is.null(snps)){
    allsnps = read_vcf(snps)
    allsnps = allsnps[lengths(allsnps$ALT)==1] ## remove multiallelic
    allsnps$ALT = unstrsplit(allsnps$ALT)
    allsnps = allsnps[nchar(allsnps$REF)==1 & nchar(allsnps$ALT) == 1] ## keep only SNPs
    othersnps = allsnps[!(allsnps %^% snps)][, c('REF', 'ALT')]
    othersnps$A = othersnps$B = othersnps$REF
    othersnps$het = FALSE
    snps = grbind(snps, othersnps)
    message('Loading and appended reference SNPs')
  }
  
  ## flag hom ALT and hom REF snps
  snps$homALT = snps$A == snps$B & snps$A == snps$ALT
  snps$homREF = snps$A == snps$B & snps$A == snps$REF
  
  ## convert snps into diploid i.e. parental haplotype coordinate
  snps.A = copy(snps); snps.B = copy(snps)
  seqlevels(snps.A) = paste(seqlevels(snps.A), "A")
  seqlevels(snps.B) = paste(seqlevels(snps.B), "B")
  values(snps.A)[, "allele"] = copy(values(snps.A)[, "A"])
  values(snps.B)[, "allele"] = copy(values(snps.B)[, "B"])
  hapsnps = grbind(snps.A, snps.B)
}



#' {
#' 
#'   
#'   parseobj = OptionParser(option_list=option_list)
#'   opt = parse_args(parseobj)
#'   
#'   if (is.null(opt$libdir) | (is.null(opt$junctions)) | is.null(opt$vcf) | is.null(opt$bias)| is.null(opt$nbias))
#'     stop(print_help(parseobj))
#'   
#'   print(opt)
#'   
#'   print(.libPaths())
#'   options(error=function() { traceback(2); quit("no", 1) })
#'   
#'   ## keep record of run
#'   writeLines(paste(paste('--', names(opt), ' ', sapply(opt, function(x) paste(x, collapse = ',')), sep = '', collapse = ' '), sep = ''), paste(opt$outdir, 'cmd.args', sep = '/'))
#'   saveRDS(opt, paste(opt$outdir, 'cmd.args.rds', sep = '/'))
#' }
#' 
#' suppressWarnings(expr = {
#'   suppressMessages(expr = {
#'     suppressPackageStartupMessages(expr = {
#'       library(skitools)
#'       library(skidb)
#'       library(data.table)
#'       library(gGnome)
#'       library(gUtils)
#'       library(bamUtils)
#'       ## devtools::load_all('~/git/gGnome')
#'     })
#'   })
#' })
#' 
#' setDTthreads(1)
#' 
#' system(paste('mkdir -p',  opt$outdir))
#' 
#' if (!is.null(opt$unmappable))
#'   CNun = opt$unmappable %>% readRDS
#' 
#' message('Loading phased SNPs')
#' snps  = read_vcf(opt$vcf, geno = TRUE)
#' fn = rev(names(values(snps)))[1]
#' snps = snps[lengths(snps$ALT)==1] ## remove multiallelic
#' snps$ALT = unstrsplit(snps$ALT)
#' snps = snps[nchar(snps$REF)==1 & nchar(snps$ALT) == 1]
#' snps = snps[, c('REF', 'ALT', fn)] %>% as.data.frame %>% .separate(fn, c('A', 'B'), '\\|') %>% dt2gr %>% gr.sub
#' if (grepl("rds$", opt$junctions)) {
#'   junctions = jJ(readRDS(opt$junctions))
#' } else {
#'   junctions = jJ(opt$junctions)
#' }
#' snps$het = snps$A != snps$B
#' snpmat = values(snps)[, c('REF', 'ALT')] %>% as.matrix
#' snps$A = snpmat[cbind(1:length(snps), as.numeric(snps$A)+1)] ## convert integer phases to bases
#' snps$B = snpmat[cbind(1:length(snps), as.numeric(snps$B)+1)]
#' message('ingested junctions and phased SNPs')
#' 
#' ## snpmat = values(snps)[, c('REF', 'ALT')] %>% as.matrix
#' ## snps$A = snpmat[cbind(1:length(snps), as.numeric(snps$A)+1)]
#' ## snps$B = snpmat[cbind(1:length(snps), as.numeric(snps$B)+1)]
#' 
#' ## additional snps that are homozogyous REF in the reference sample may not be
#' ## present in the phased VCF, so we can pull these from a reference DB e.g. hapmap
#' if (!is.null(opt$snps))
#' {
#'   allsnps = read_vcf(opt$snps)
#'   allsnps = allsnps[lengths(allsnps$ALT)==1] ## remove multiallelic
#'   allsnps$ALT = unstrsplit(allsnps$ALT)
#'   allsnps = allsnps[nchar(allsnps$REF)==1 & nchar(allsnps$ALT) == 1] ## keep only SNPs
#'   othersnps = allsnps[!(allsnps %^% snps)][, c('REF', 'ALT')]
#'   othersnps$A = othersnps$B = othersnps$REF
#'   othersnps$het = FALSE
#'   snps = grbind(snps, othersnps)
#'   message('Loading and appended reference SNPs')
#' }
#' 
#' ## flag hom ALT and hom REF snps
#' snps$homALT = snps$A == snps$B & snps$A == snps$ALT
#' snps$homREF = snps$A == snps$B & snps$A == snps$REF
#' 
#' ## convert snps into diploid i.e. parental haplotype coordinate
#' snps.A = copy(snps); snps.B = copy(snps)
#' seqlevels(snps.A) = paste(seqlevels(snps.A), "A")
#' seqlevels(snps.B) = paste(seqlevels(snps.B), "B")
#' values(snps.A)[, "allele"] = copy(values(snps.A)[, "A"])
#' values(snps.B)[, "allele"] = copy(values(snps.B)[, "B"])
#' hapsnps = grbind(snps.A, snps.B)
#' ## hapsnps = grbind(
#' ##   gr2dt(snps)[, seqnames := paste(seqnames, 'A')][, allele := A] %>% dt2gr %>% cc('allele'),
#' ##   gr2dt(snps)[, seqnames := paste(seqnames, 'B')][, allele := B] %>% dt2gr %>% cc('allele')
#' ## )
#' 
#' ## make new haplotype seqlengths i.e. across A and B haplotypes
#' sl = expand.grid(sl = seqlevels(snps), hap = c('A', 'B')) %>% as.data.table %>% cc(structure(seqlengths(snps)[sl], names = paste(sl, hap)))
#' 
#' jjhap = junctions
#' if (length(junctions))
#' {
#'   ## initial graph to make eclusters
#'   igg = gG(junctions = junctions)
#'   igg$eclusters(thresh = 1e5)
#'   ## only ALT junctions
#'   jj = igg$junctions[type == "ALT"][, 'ecluster']
#'   jj = jj$set(ecluster = ifelse(is.na(jj$dt$ecluster), max(c(0, jj$dt$ecluster), na.rm = TRUE) + 1:length(jj), jj$dt$ecluster))
#'   
#'   ## now assign every ecluster and chromosome the same haplotype and jcn given some ploidy
#'   hapdt = (jj$grl %>% grl.unlist) %>% gr2dt %>% cc(.(seqnames, ecluster)) %>% unique %>% cc(hap := sample(c('A', 'B'), .N, replace = TRUE), by = .(seqnames, ecluster)) %>% setkeyv(c('seqnames', 'ecluster'))
#'   
#'   hapdt$jcn = rexp(nrow(hapdt), rate = 2/opt$tau) %>% ceiling
#'   jj$set(hap1 = hapdt[.(seqnames(jj$left) %>% as.character, jj$dt$ecluster), hap])
#'   jj$set(hap2 = hapdt[.(seqnames(jj$right) %>% as.character, jj$dt$ecluster), hap])
#'   jj$set(cn = hapdt[.(seqnames(jj$right) %>% as.character, jj$dt$ecluster), jcn])
#'   
#'   ## now lift junctions into haplotype coordinates
#'   left = jj$left %>% gr2dt %>% cc(seqnames := paste(seqnames, jj$dt$hap1)) %>% dt2gr(seqlengths = sl)
#'   right = jj$right %>% gr2dt %>% cc(seqnames := paste(seqnames, jj$dt$hap2)) %>% dt2gr(seqlengths = sl)
#'   hapgrl = split(grbind(left, right), rep(1:length(left), 2))
#'   values(hapgrl) = jj$dt
#'   jjhap = jJ(hapgrl)
#'   message('assigned junctions to germline haplotypes')
#'   
#'   ## make genome graph with lb 1 on all provided junctions
#'   jjhap$set(cn = NA)
#'   jjhap$set(lb = 1)
#' }
#' 
#' 
#' breaks = NULL
#' if (!is.null(opt$unmappable))
#' {
#'   CNun = rbind(gr2dt(CNun)[, seqnames := paste(seqnames, 'A')],
#'                gr2dt(CNun)[, seqnames := paste(seqnames, 'B')]) %>% dt2gr
#'   breaks = gr.sample(CNun, opt$numbreaks) ## sample some random unmappable breaks
#' }
#' 
#' #' zchoo Monday, Apr 25, 2022 11:19:38 AM
#' ## add option for including a CNLOH junction
#' if (opt$cnloh) {
#'   
#'   message("Creating a CNLOH edge")
#'   
#'   ## get the gaps between current breakends
#'   ## in **diploid** coordinates
#'   if (grepl("rds$", opt$junctions)) {
#'     og.junctions = jJ(readRDS(opt$junctions))
#'   } else {
#'     og.junctions = jJ(opt$junctions)
#'   }
#'   
#'   if (length(og.junctions)) {
#'     jjhap.bnds = unlist(og.junctions$grl)
#'   } else {
#'     jjhap.bnds = GRanges(seqinfo = seqinfo(og.junctions), seqlengths = seqlengths(og.junctions))
#'   }
#'   
#'   if (is.null(breaks) || is.na(breaks) || length(breaks) == 0) {
#'     bnds = gr.stripstrand(jjhap.bnds)
#'   } else {
#'     bnds = grbind(gr.stripstrand(jjhap.bnds), gr.stripstrand(breaks))
#'   }
#'   
#'   ## get gaps that are at least 1e5 from the breaks
#'   gaps = gaps(bnds + 1e5)
#'   gaps = gaps %Q% (strand(gaps) == "*")
#'   ## only autosomes
#'   gaps = gaps %Q% (as.character(seqnames(gaps)) %in% as.character(1:22))
#'   
#'   ## sample a point (one for now but make this adjustable in the future?)
#'   ## this would need to be iterative to avoid sampling points taht are too close together
#'   cnloh.br = gr.sample(gaps, k = 1, wid = 1)
#'   
#'   message("Sampled a CNLOH edge located at: ",
#'           paste(seqnames(cnloh.br), ":", GenomicRanges::start(cnloh.br)))
#'   
#'   ## create a junction from haplotype A to B
#'   ## but should in theory make a vector giving the A/B parity
#'   cnloh.br.haploid = GRangesList(GRanges(seqnames = c(paste(seqnames(cnloh.br), "A"),
#'                                                       paste(seqnames(cnloh.br), "B")),
#'                                          ranges = IRanges(start = c(GenomicRanges::start(cnloh.br),
#'                                                                     GenomicRanges::start(cnloh.br) + 1),
#'                                                           width = 1),
#'                                          strand = c("-", "+")))
#'   
#'   jjcnloh = jJ(cnloh.br.haploid)
#'   jjcnloh$set(cn = NA)
#'   jjcnloh$set(lb = 1) ## force in this junction
#'   jjcnloh$set(ub = 1) ## force in this junction
#'   ## mark this in some way?
#'   jjcnloh$set(cnloh = TRUE)
#'   jjhap = c(jjhap, jjcnloh)
#' }
#' 
#' gg = gG(junctions = jjhap, breaks = breaks)
#' gg$nodes$mark(cn = opt$tau/2)
#' 
#' ## if using cnloh we will need to fix some nodes
#' ## hmm let's fix the start to be equal?
#' if (opt$cnloh) {
#'   message("Fixing nodes to be equal...")
#'   n1 = gg$nodes$dt[seqnames == paste(seqnames(cnloh.br), "A") & end == GenomicRanges::start(cnloh.br), node.id]
#'   n2 = gg$nodes$dt[seqnames == paste(seqnames(cnloh.br), "B") & end == GenomicRanges::start(cnloh.br), node.id]
#'   gg$nodes[node.id == n1]$mark(ub = ceiling(opt$tau/2), lb = ceiling(opt$tau/2))
#'   gg$nodes[node.id == n2]$mark(ub = ceiling(opt$tau/2), lb = ceiling(opt$tau/2))
#' }
#' 
#' gg$nodes$mark(weight = 0.1*gg$nodes$dt$width/1e6)
#' message('built basic graph from junctions and copy number breaks')
#' 
#' 
#' ## fitted haplotype graph, where every provided alt edge is given some nonzero copy number
#' saveRDS(gg, paste0(opt$outdir, '/tmp.gg.rds'))
#' ggb = balance(gg, lambda = 10000, epgap = 1e-6, ism = FALSE, verbose = 2, tilim = 600)
#' ggb$nodes$mark(hap = gg$nodes$dt$seqnames %>% as.character %>% strsplit(' ') %>% sapply('[', 2))
#' message('balanced phased graph constraining all junctions to be incorporated yielding diploid ploidy ', gGnome:::ploidy(ggb) %>% signif(3))
#' 
#' ## lift haplotype graph from diploid coordinates (i.e. separate coordiante for each haplotype)
#' ## back to standard haploid genome coordinates
#' ## there should be two nodes per haploid coordinate here
#' ## lift phased graph to haploid coordinates
#' ## use just the final space?
#' nodes = ggb$nodes$gr %>% gr2dt %>% .separate(old = 'seqnames', sep = ' (?=[^ ]+$)', new = c('seqnames', 'hap'), perl = TRUE) %>% dt2gr
#' ## ggl = gG(nodes = nodes, edges = ggb$edges$dt) %>% gGnome:::loosefix ## pipe and ":::" not playing nice
#' ggl = gGnome:::loosefix(gG(nodes = nodes, edges = ggb$edges$dt))
#' ggl$edges$mark(type = ifelse(ggl$edges$class == 'REF', 'REF', 'ALT'))
#' ggl$set(y.field = 'cn', purity = opt$alpha, ploidy = gGnome:::ploidy(ggl))
#' message('lifted phased graph from diploid to haploid coordinates')
#' 
#' ## also make disjoined / collapsed jabba graph, will be one of the key outputs
#' ## there will be one node per haploid coordinate here
#' ggd = ggl$copy
#' ggd$nodes$mark(acn = ifelse(ggd$nodes$dt$hap == 'A', ggd$nodes$dt$cn, 0))
#' ggd$nodes$mark(bcn = ifelse(ggd$nodes$dt$hap == 'B', ggd$nodes$dt$cn, 0))
#' ggd = ggd$disjoin()
#' ggd$edges$mark(type = ifelse(ggd$edges$class == 'REF', 'REF', 'ALT'))
#' ggd$set(y.field = 'cn', purity = opt$alpha, ploidy = gGnome:::ploidy(ggd))
#' message('made collapsed jabba style graph with one node per haploid coordinate, keeping track of allelic copy numbers giving ploidy ', gGnome:::ploidy(ggd) %>% signif(3))
#' 
#' ## now sample from ggb and ggd to get total and het coverage
#' bias = opt$bias %>% readRDS %>% gr.nochr
#' fn = names(values(bias))[1] ## take first column as value
#' #bias = bias %>% rebin(field = fn, 1e4) ## smooth out across 10kb temporarily silenced, addy
#' bias$bias = values(bias)[[1]]/mean(values(bias)[[1]], na.rm = TRUE)
#' message('ingested and processed tumor bias GRanges')
#' 
#' ## process second normal sample which will be used to simulate the normal depth
#' nbias = opt$nbias %>% readRDS %>% gr.nochr
#' fn = names(values(nbias))[1]
#' #nbias = nbias %>% rebin(field = fn, 1e4)#temporarily silenced, addy
#' nbias$bias = values(nbias)[[1]]/mean(values(nbias)[[1]], na.rm = TRUE)
#' message('ingested and processed normal bias GRanges')
#' 
#' cov = .simcov(ggd$nodes$gr, purity = opt$alpha, basecov = opt$coverage, binsize = opt$width, normalize = FALSE, bias = bias[, 'bias'], poisson = opt$poisson)
#' tmpgr = ggd$nodes$gr; tmpgr$cn = 2
#' ncov = .simcov(tmpgr, purity = 1, basecov = opt$ncoverage, normalize = FALSE, binsize = opt$width, bias = nbias[, 'bias'], poisson = opt$poisson)
#' 
#' ## final output binned coverage
#' outcov = cov
#' outcov$tumor = cov$cov
#' outcov$normal = ncov$cov
#' outcov$ratio = outcov$tumor/outcov$normal
#' outcov$ratio = outcov$ratio/median(outcov$ratio, na.rm = TRUE)
#' message('simulated binned total coverage for tumor and matched normal')
#' 
#' ## simulate SNP coverage in diploid genome
#' 
#' ## lift bias GRanges to diploid coordinates
#' a.bias = copy(bias)
#' b.bias = copy(bias)
#' seqlevels(a.bias) = paste(seqlevels(a.bias), "A")
#' seqlevels(b.bias) = paste(seqlevels(b.bias), "B")
#' hapbias = grbind(a.bias, b.bias)
#' ## hapbias = rbind(as.data.table(bias)[, seqnames := paste(seqnames, 'A')],
#' ##                 as.data.table(bias)[, seqnames := paste(seqnames, 'B')]) %>% dt2gr
#' 
#' hapsnpcov = .simcov(ggb$nodes$gr, bins = hapsnps, bias = hapbias, diploid = FALSE, basecov = opt$coverage/2, normalize = FALSE)
#' hapsnpcov$count = round(hapsnpcov$cov)
#' 
#' ## now populate snps
#' snpcov = hapsnpcov[, c("count", "allele", "cn")] %>% as.data.table %>% .separate('seqnames', c('seqnames', 'hap'), ' ') %>%  dcast.data.table(seqnames + start + end ~ hap, value.var = c('count', 'cn', 'allele')) %>% dt2gr
#' 
#' snpcov$REF = snps$REF[gr.match(snpcov, snps)]
#' snpcov$ALT = snps$ALT[gr.match(snpcov, snps)]
#' snpcov$count_ALT = snpcov$count_A*sign(snpcov$allele_A == snpcov$ALT) + snpcov$count_B*sign(snpcov$allele_B == snpcov$ALT)
#' snpcov$count_REF = snpcov$count_A*sign(snpcov$allele_A == snpcov$REF) + snpcov$count_B*sign(snpcov$allele_B == snpcov$REF)
#' snpcov$cn_ALT = snpcov$cn_A*sign(snpcov$allele_A == snpcov$ALT) + snpcov$cn_B*sign(snpcov$allele_B == snpcov$ALT)
#' snpcov$cn_REF = snpcov$cn_A*sign(snpcov$allele_A == snpcov$REF) + snpcov$cn_B*sign(snpcov$allele_B == snpcov$REF)
#' snpcov$high_count = pmax(snpcov$count_ALT, snpcov$count_REF)
#' snpcov$low_count = pmin(snpcov$count_ALT, snpcov$count_REF)
#' snpcov$het = snps$het[gr.match(snpcov, snps)]
#' snpcov$homALT = snps$homALT[gr.match(snpcov, snps)]
#' snpcov$homREF = snps$homREF[gr.match(snpcov, snps)]
#' message('simulated phased tumor allelic coverage and lifted to haploid coordinates')
#' 
#' ## generate normal sample hets with normal bias
#' snpcov = snpcov %$% nbias[, 'bias']
#' snpcov$ncount_ALT = rpois(length(snpcov), (snpcov$homALT + snpcov$het) * opt$coverage/2 * snpcov$bias)
#' snpcov$ncount_REF = rpois(length(snpcov), (snpcov$homREF + snpcov$het) * opt$ncoverage/2 * snpcov$bias)
#' message('simulated normal allelic coverage')
#' 
#' ## prepare het_pileups_wgs - like outtput
#' hets = as.data.table(snpcov)[, .(seqnames, start, end, strand,
#'                                  Tumor_Seq_Allele1 = ALT,
#'                                  Reference_Allele = REF,
#'                                  alt.count.t = count_ALT,
#'                                  ref.count.t = count_REF,
#'                                  alt.frac.t = ifelse(count_ALT > 0,
#'                                                      count_ALT / (count_ALT + count_REF),
#'                                                      0),
#'                                  ref.frac.t = ifelse(count_REF > 0,
#'                                                      count_REF / (count_ALT + count_REF),
#'                                                      0),
#'                                  alt.count.n = as.numeric(ncount_ALT),
#'                                  ref.count.n = as.numeric(ncount_REF))]
#' 
#' 
#' hets[, ref.frac.n := ifelse(ref.count.n > 0, ref.count.n / (ref.count.n + alt.count.n), 0)]
#' hets[, alt.frac.n := ifelse(alt.count.n > 0, alt.count.n / (ref.count.n + alt.count.n), 0)]
#' 
#' ## make gTracks
#' redblue = rbind(gr2dt(snpcov)[het == TRUE, count := high_count][, col := alpha('red', 0.1)],
#'                 gr2dt(snpcov)[het == TRUE, count := low_count][, col := alpha('blue', 0.1)]) %>% dt2gr
#' 
#' purplegreen = rbind(gr2dt(snpcov)[het == TRUE, count := count_A][, col := alpha('purple', 0.1)],
#'                     gr2dt(snpcov)[het == TRUE, count := count_B][, col := alpha('green', 0.1)]) %>% dt2gr
#' 
#' 
#' gt.rb = gTrack(redblue, y.field = 'count', y0 = 0, circle = TRUE, y1 = opt$coverage*2, lwd.border = 0.2, name = 'high-low', max.ranges = 1e4)
#' gt.pg = gTrack(purplegreen, y.field = 'count', y0 = 0, circle = TRUE, y1 = opt$coverage*2, lwd.border = 0.2, name = 'phased', max.ranges = 1e4)
#' 
#' #' zchoo Monday, Apr 25, 2022 12:44:10 PM
#' ## restrict the max.ranges on coverage
#' gt.cov = gTrack(outcov, y.field = c('tumor'), lwd.border = 0.2, circle = TRUE, max.ranges = 1e4)
#' gt.ncov = gTrack(outcov, y.field = c('normal'), lwd.border = 0.2, circle = TRUE, max.ranges = 1e4)
#' gt.tcov = gTrack(outcov, y.field = c('ratio'), lwd.border = 0.2, circle = TRUE, y1 = 3, y0 = 0, max.ranges = 1e4)
#' tmp.gg = ggl$copy
#' tmp.gg$nodes$mark(cn = ifelse(tmp.gg$nodes$gr$hap == 'A', tmp.gg$nodes$gr$cn + 0.2, tmp.gg$nodes$gr$cn - 0.2))
#' tmp.gg$nodes$mark(col = ifelse(tmp.gg$nodes$gr$hap == 'A', 'purple', 'green'))
#' gt = c(gt.rb, gt.pg, gt.ncov, gt.cov, gt.tcov, tmp.gg$gtrack(name = 'phased'), ggd$gtrack(name = 'collapsed'))
#' message('made gTracks')
#' 
#' fwrite(hets, paste0(opt$outdir, '/sites.txt'))
#' saveRDS(outcov, paste0(opt$outdir, '/cov.rds'))
#' saveRDS(ggl, paste0(opt$outdir, '/graph.phased.rds'))
#' saveRDS(ggd, paste0(opt$outdir, '/graph.unphased.rds'))
#' saveRDS(snpcov, paste0(opt$outdir, '/snps.rds'))
#' saveRDS(gt, paste0(opt$outdir, '/gt.rds'))
#' message('dumped out files and finished')
#' }, echo = FALSE)