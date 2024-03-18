#' @import skitools
#' @import data.table
#' @import gGnome
#' @import gUtils
#' @import bamUtils
#' @import skidb
#' @import gTrack
#' @import dplyr
#' @importMethodsFrom GenomicRanges as.data.frame

#' @name ggsim
#' 
#' @description function to simulate a gGraph complete with SNPs and junction phasing, and coverage
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
#' @param tau Target ploidy
#' @param poisson add poisson noise
#' @param numbreaks number of additional breaks to provide in cn unmappable regions
#' @param width Binwidth for binned read depth
#' @param cnloh add a cnloh edge
#' @param standard.chr chromosomes to simulate
#' @param par.path path to pseudoautosomal regions on X, Y chromosome
#' @param outdir output directory
#'
#' @export
#'
#' @examples
ggsim = function(junctions,
                 vcf,
                 bias,
                 nbias,
                 snp = NULL,
                 unmappable = NULL,
                 coverage,
                 ncoverage,
                 alpha,
                 tau,
                 poisson = TRUE,
                 numbreaks,
                 width = 1000,
                 cnloh = FALSE,
                 standard.chr = c(1:22, "X", "Y"),
                 outdir,
                 par.path = system.file("extdata", "PAR_hg19.rds", package = 'ggSim'))
{
  if (is.null(junctions) | is.null(vcf) | is.null(bias)| is.null(nbias))
    stop("Check junctions, vcf, bias, nbias to ensure they are not empty inputs.")
  
  setDTthreads(1)
  
  system(paste('mkdir -p',  outdir))

  if (!is.null(unmappable))
    CNun = unmappable %>% readRDS

  par.file = par.path %>% readRDS()
  
  bias = bias %>% readRDS %>% gr.nochr
  fn = names(values(bias))[1] ## take first column as value
  #bias = bias %>% rebin(field = fn, 1e4) ## smooth out across 10kb temporarily silenced, addy
  bias$bias = values(bias)[[1]]/mean(values(bias)[[1]], na.rm = TRUE)
  bias$autosome = (bias %^% par.file) | (gr2dt(bias)$seqnames %in% c(1:22)) #tag (pseudo/)autosomal
  message('ingested and processed tumor bias GRanges')
  
  ## process second normal sample which will be used to simulate the normal depth
  nbias = nbias %>% readRDS %>% gr.nochr
  fn = names(values(nbias))[1]
  #nbias = nbias %>% rebin(field = fn, 1e4)#temporarily silenced, addy
  nbias$bias = values(nbias)[[1]]/mean(values(nbias)[[1]], na.rm = TRUE)
  nbias$autosome = (nbias %^% par.file) | (gr2dt(nbias)$seqnames %in% c(1:22)) #tag (pseudo/)autosomal
  message('ingested and processed normal bias GRanges')
  
  #make the haploid coverages diploid on the X and Y chromosome... this will matter for sim coverages later
  bias = gr2dt(bias)[seqnames %in% c("X", "Y"), bias := ifelse(bias == 0, 0, bias / mean(bias)), by = .(autosome, seqnames)] %>% dt2gr
  nbias = gr2dt(nbias)[seqnames %in% c("X", "Y"), bias := ifelse(bias == 0, 0, bias / mean(bias)), by = .(autosome, seqnames)] %>% dt2gr
  
  message('Loading phased SNPs')
  snps  = skidb::read_vcf(vcf, geno = TRUE, verbose = TRUE)
  fn = rev(names(values(snps)))[1]
  snps = snps[lengths(snps$ALT)==1] ## remove multiallelic
  snps$ALT = unstrsplit(snps$ALT)
  snps = snps %Q% (Biostrings::nchar(REF) == 1 && Biostrings::nchar(ALT) == 1)
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
  snps = snps %Q% (seqnames %in% standard.chr)
  seqlevels(snps) <- seqlevels(snps)[seqlevels(snps) %in% standard.chr]
  message('ingested junctions and phased SNPs')
  
  ## sex typing based on input vcf
  sex <- c("M", "F")[any(snps[seqnames(snps) == "X"]$het) %>% as.numeric() + 1]
  message(sprintf("Your input vcf was determined to be %s. Sex of sim genome will be %s.", sex, sex))
  
  #additional snps that are homozogyous REF in the reference sample may not be
  #present in the phased VCF, so we can pull these from a reference DB e.g. hapmap
  if (!is.null(snp)){
    allsnps = read_vcf(snp)
    allsnps = allsnps[lengths(allsnps$ALT)==1] ## remove multiallelic
    allsnps$ALT = unstrsplit(allsnps$ALT)
    allsnps = allsnps[Biostrings::nchar(allsnps$REF)==1 & Biostrings::nchar(allsnps$ALT) == 1] ## keep only SNPs
    othersnps = allsnps[!(allsnps %^% snps)][, c('REF', 'ALT')]
    othersnps$A = othersnps$B = othersnps$REF
    othersnps$het = FALSE
    othersnps@ranges@NAMES <- NULL
    othersnps = othersnps %Q% (seqnames %in% standard.chr)
    seqlevels(othersnps) <- seqlevels(othersnps)[seqlevels(othersnps) %in% standard.chr]
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
  
  #get rid of alleles based on sex
  hapsnps = grbind(snps.A, snps.B)
  if(sex == "M"){
    hapsnps = hapsnps %Q% (!seqnames %in% c("X B", "Y B"))
    seqlevels(hapsnps) = seqlevels(hapsnps)[!seqlevels(hapsnps) %in% c("X B", "Y B")] } else { hapsnps = hapsnps %Q% (!seqnames %in% c("Y A", "Y B"))
    seqlevels(hapsnps) = seqlevels(hapsnps)[!seqlevels(hapsnps) %in% c("Y A", "Y B")]}
      
  ## make new haplotype seqlengths i.e. across A and B haplotypes
  sl = seqlengths(hapsnps)
  
  jjhap = junctions
  if (length(junctions)) {
    ## initial graph to make eclusters
    igg = gG(junctions = junctions)
    igg$eclusters(thresh = 1e5)
    
    ## only ALT junctions
    jj = igg$junctions[type == "ALT"][, 'ecluster']
    #give NA junctions a unique cluster
    jj = jj$set(ecluster = ifelse(is.na(jj$dt$ecluster), max(c(0, jj$dt$ecluster), na.rm = TRUE) + 1:length(jj), jj$dt$ecluster)) 

    chrom_allele = hapsnps %>% gr2dt %>% cc(.(seqnames)) %>% unique %>% transmute(chrom = gsub(" .*", "", seqnames), allele = gsub(".* ", "", seqnames))
    
    ## now assign every ecluster and chromosome the same haplotype and jcn given some ploidy
    hapdt = (jj$grl %>% grl.unlist) %>% 
      gr2dt %>% 
      cc(.(seqnames, ecluster)) %>% 
      unique %>%
      cc(hap := sample(c('A', 'B'), .N, replace = TRUE), by = .(seqnames, ecluster)) %>%
      rowwise() %>%
      mutate(hap = ifelse(!hap %in% chrom_allele[chrom == seqnames]$allele, "A", hap))%>%
      data.table() %>%
      setkeyv(c('seqnames', 'ecluster'))

    hapdt$jcn = rexp(nrow(hapdt), rate = 2/tau) %>% ceiling
    jj$set(hap1 = hapdt[.(seqnames(jj$left) %>% as.character, jj$dt$ecluster), hap])
    jj$set(hap2 = hapdt[.(seqnames(jj$right) %>% as.character, jj$dt$ecluster), hap])
    jj$set(cn = hapdt[.(seqnames(jj$right) %>% as.character, jj$dt$ecluster), jcn])
    
    
    ## now lift junctions into haplotype coordinates
    left = jj$left %>% gr2dt 
    left$hap = jj$dt$hap1
    left = left %>% cc(seqnames := paste(seqnames, hap)) %>% dt2gr(seqlengths = sl)
    #%>% mutate(seqnames = paste(seqnames, left.hap)) %>% dt2gr(seqlengths = sl)
    right = jj$right %>% gr2dt 
    right$hap = jj$dt$hap2
    right = right %>% cc(seqnames := paste(seqnames, hap)) %>% dt2gr(seqlengths = sl)
    hapgrl = GenomicRanges::split(grbind(left, right), rep(1:length(left), 2))
    #hapgrl@metadata <- jj$dt
    values(hapgrl) <- jj$dt
    jjhap = jJ(hapgrl)
    message('assigned junctions to germline haplotypes')

    ## make genome graph with lb 1 on all provided junctions
    jjhap$set(cn = NA)
    jjhap$set(lb = 1)
  }
  
  breaks = NULL
  if (!is.null(unmappable)) {
    CNun = unmappable %>% readRDS
    CNun = rbind(gr2dt(CNun)[, seqnames := paste(seqnames, 'A')],
                gr2dt(CNun)[, seqnames := paste(seqnames, 'B')]) %>% 
      dt2gr %Q% (seqnames %in% names(sl))
    seqlevels(CNun) <- names(sl)
    breaks = gr.sample(CNun, numbreaks)
  }
  
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
 
  gg = gG(junctions = jjhap, breaks = breaks)
  gg$nodes$mark(cn = tau/2)
  
  #' ## if using cnloh we will need to fix some nodes
  #' ## hmm let's fix the start to be equal?
  #' if (cnloh) {
  #'   message("Fixing nodes to be equal...")
  #'   n1 = gg$nodes$dt[seqnames == paste(seqnames(cnloh.br), "A") & end == GenomicRanges::start(cnloh.br), node.id]
  #'   n2 = gg$nodes$dt[seqnames == paste(seqnames(cnloh.br), "B") & end == GenomicRanges::start(cnloh.br), node.id]
  #'   gg$nodes[node.id == n1]$mark(ub = ceiling(opt$tau/2), lb = ceiling(opt$tau/2))
  #'   gg$nodes[node.id == n2]$mark(ub = ceiling(opt$tau/2), lb = ceiling(opt$tau/2))
  #' }
  
  gg$nodes$mark(weight = gg$nodes$dt$width/1e7)
  message('built basic graph from junctions and copy number breaks')

  #' ## fitted haplotype graph, where every provided alt edge is given some nonzero copy number
  #' saveRDS(gg, paste0(outdir, '/tmp.gg.rds'))
  ggb = balance(gg, 
                lambda = 10000, 
                epgap = 1e-6, 
                ism = FALSE, 
                verbose = 2, 
                tilim = 600)
  ggb$nodes$mark(hap = gg$nodes$dt$seqnames %>% as.character %>% strsplit(' ') %>% sapply('[', 2))
  message('balanced phased graph constraining all junctions to be incorporated yielding diploid ploidy ', gGnome:::ploidy(ggb) %>% signif(3))
 
  ## lift haplotype graph from diploid coordinates (i.e. separate coordiante for each haplotype)
  ## back to standard haploid genome coordinates
  ## there should be two nodes per haploid coordinate here
  ## lift phased graph to haploid coordinates
  ## use just the final space?
  nodes = ggb$nodes$gr %>% gr2dt %>% .separate(old = 'seqnames', 
                                               sep = ' (?=[^ ]+$)', 
                                               new = c('seqnames', 'hap'), 
                                               perl = TRUE) %>%
    dt2gr
  ## ggl = gG(nodes = nodes, edges = ggb$edges$dt) %>% gGnome:::loosefix ## pipe and ":::" not playing nice
  ggl = gGnome:::loosefix(gG(nodes = nodes, edges = ggb$edges$dt))
  ggl$edges$mark(type = ifelse(ggl$edges$class == 'REF', 'REF', 'ALT'))
  ggl$set(y.field = 'cn', purity = alpha, ploidy = gGnome:::ploidy(ggl))
  message('lifted phased graph from diploid to haploid coordinates')
  
  ## also make disjoined / collapsed jabba graph, will be one of the key outputs
  ## there will be one node per haploid coordinate here
  ggd = ggl$copy
  ggd$nodes$mark(acn = ifelse(ggd$nodes$dt$hap == 'A', ggd$nodes$dt$cn, 0))
  ggd$nodes$mark(bcn = ifelse(ggd$nodes$dt$hap == 'B', ggd$nodes$dt$cn, 0))
  ggd = ggd$disjoin()
  ggd$edges$mark(type = ifelse(ggd$edges$class == 'REF', 'REF', 'ALT'))
  ggd$set(y.field = 'cn', purity = alpha, ploidy = gGnome:::ploidy(ggd))
  message('made collapsed jabba style graph with one node per haploid coordinate, keeping track of allelic copy numbers giving ploidy ', gGnome:::ploidy(ggd) %>% signif(3))

  ###########################################################################
  
  ## now sample from ggb and ggd to get total and het coverage

  cov = simulate_coverage(gr = ggd$nodes$gr,
                          purity = alpha,
                          basecov = coverage,
                          binsize = width,
                          normalize = F,
                          bias = bias[,'bias'],
                          poisson = poisson)
  tmpgr = ggd$nodes$gr; tmpgr$cn = 2
  ncov = simulate_coverage(gr = ggd$nodes$gr, 
                           purity = 1,
                           basecov = ncoverage,
                           normalize = FALSE,
                           binsize = width,
                           bias = nbias[, 'bias'],
                           poisson = poisson)
  
  ## final output binned coverage
  outcov = cov
  outcov$tumor = cov$cov
  outcov$normal = ncov$cov
  outcov$ratio = outcov$tumor/outcov$normal
  outcov$ratio = outcov$ratio/median(outcov$ratio, na.rm = TRUE)
  message('simulated binned total coverage for tumor and matched normal')

  ## simulate SNP coverage in diploid genome

  ## lift bias GRanges to diploid coordinates
  a.bias = copy(bias)
  b.bias = copy(bias)
  seqlevels(a.bias) = paste(seqlevels(a.bias), "A")
  seqlevels(b.bias) = paste(seqlevels(b.bias), "B")
  hapbias = grbind(a.bias, b.bias)
  ## hapbias = rbind(as.data.table(bias)[, seqnames := paste(seqnames, 'A')],
  ##                 as.data.table(bias)[, seqnames := paste(seqnames, 'B')]) %>% dt2gr

  hapsnpcov = simulate_coverage(gr = ggb$nodes$gr,
                                bins = hapsnps,                      
                                bias = hapbias,
                                diploid = FALSE,
                                basecov = coverage/2,
                                normalize = FALSE)
  hapsnpcov$count = round(hapsnpcov$cov)

  ## now populate snps
  snpcov = hapsnpcov[, c("count", "allele", "cn")] %>% as.data.table %>% .separate('seqnames', c('seqnames', 'hap'), ' ') %>%  dcast.data.table(seqnames + start + end ~ hap, value.var = c('count', 'cn', 'allele')) %>% dt2gr

  snpcov$REF = snps$REF[gr.match(snpcov, snps)]
  snpcov$ALT = snps$ALT[gr.match(snpcov, snps)]
  snpcov$count_ALT = snpcov$count_A*sign(snpcov$allele_A == snpcov$ALT) + snpcov$count_B*sign(snpcov$allele_B == snpcov$ALT)
  snpcov$count_REF = snpcov$count_A*sign(snpcov$allele_A == snpcov$REF) + snpcov$count_B*sign(snpcov$allele_B == snpcov$REF)
  snpcov$cn_ALT = snpcov$cn_A*sign(snpcov$allele_A == snpcov$ALT) + snpcov$cn_B*sign(snpcov$allele_B == snpcov$ALT)
  snpcov$cn_REF = snpcov$cn_A*sign(snpcov$allele_A == snpcov$REF) + snpcov$cn_B*sign(snpcov$allele_B == snpcov$REF)
  snpcov$high_count = pmax(snpcov$count_ALT, snpcov$count_REF)
  snpcov$low_count = pmin(snpcov$count_ALT, snpcov$count_REF)
  snpcov$het = snps$het[gr.match(snpcov, snps)]
  snpcov$homALT = snps$homALT[gr.match(snpcov, snps)]
  snpcov$homREF = snps$homREF[gr.match(snpcov, snps)]
  message('simulated phased tumor allelic coverage and lifted to haploid coordinates')

  ## generate normal sample hets with normal bias
  snpcov = snpcov %$% nbias[, 'bias']
  snpcov$ncount_ALT = rpois(length(snpcov), (snpcov$homALT + snpcov$het) * coverage/2 * snpcov$bias)
  snpcov$ncount_REF = rpois(length(snpcov), (snpcov$homREF + snpcov$het) * ncoverage/2 * snpcov$bias)
  message('simulated normal allelic coverage')

  ## prepare het_pileups_wgs - like outtput
  hets = as.data.table(snpcov)[, .(seqnames, start, end, strand,
                                   Tumor_Seq_Allele1 = ALT,
                                   Reference_Allele = REF,
                                   alt.count.t = count_ALT,
                                   ref.count.t = count_REF,
                                   alt.frac.t = ifelse(count_ALT > 0,
                                                       count_ALT / (count_ALT + count_REF),
                                                       0),
                                   ref.frac.t = ifelse(count_REF > 0,
                                                       count_REF / (count_ALT + count_REF),
                                                       0),
                                   alt.count.n = as.numeric(ncount_ALT),
                                   ref.count.n = as.numeric(ncount_REF))]


  hets[, ref.frac.n := ifelse(ref.count.n > 0, ref.count.n / (ref.count.n + alt.count.n), 0)]
  hets[, alt.frac.n := ifelse(alt.count.n > 0, alt.count.n / (ref.count.n + alt.count.n), 0)]

  ## make gTracks
  redblue = rbind(gr2dt(snpcov)[het == TRUE, count := high_count][, col := alpha('red', 0.1)],
                  gr2dt(snpcov)[het == TRUE, count := low_count][, col := alpha('blue', 0.1)]) %>% dt2gr

  purplegreen = rbind(gr2dt(snpcov)[het == TRUE, count := count_A][, col := alpha('purple', 0.1)],
                      gr2dt(snpcov)[het == TRUE, count := count_B][, col := alpha('green', 0.1)]) %>% dt2gr


  gt.rb = gTrack(redblue, y.field = 'count', y0 = 0, circle = TRUE, y1 = coverage*2, lwd.border = 0.2, name = 'high-low', max.ranges = 1e4)
  gt.pg = gTrack(purplegreen, y.field = 'count', y0 = 0, circle = TRUE, y1 = coverage*2, lwd.border = 0.2, name = 'phased', max.ranges = 1e4)

  #' zchoo Monday, Apr 25, 2022 12:44:10 PM
  ## restrict the max.ranges on coverage
  gt.cov = gTrack(outcov, y.field = c('tumor'), lwd.border = 0.2, circle = TRUE, max.ranges = 1e4)
  gt.ncov = gTrack(outcov, y.field = c('normal'), lwd.border = 0.2, circle = TRUE, max.ranges = 1e4)
  gt.tcov = gTrack(outcov, y.field = c('ratio'), lwd.border = 0.2, circle = TRUE, y1 = 3, y0 = 0, max.ranges = 1e4)
  tmp.gg = ggl$copy
  tmp.gg$nodes$mark(cn = ifelse(tmp.gg$nodes$gr$hap == 'A', tmp.gg$nodes$gr$cn + 0.2, tmp.gg$nodes$gr$cn - 0.2))
  tmp.gg$nodes$mark(col = ifelse(tmp.gg$nodes$gr$hap == 'A', 'purple', 'green'))
  gt = c(gt.rb, gt.pg, gt.ncov, gt.cov, gt.tcov, tmp.gg$gtrack(name = 'phased'), ggd$gtrack(name = 'collapsed'))
  message('made gTracks')

  fwrite(hets, paste0(outdir, '/sites.txt'))
  saveRDS(outcov, paste0(outdir, '/cov.rds'))
  saveRDS(ggl, paste0(outdir, '/graph.phased.rds'))
  saveRDS(ggd, paste0(outdir, '/graph.unphased.rds'))
  saveRDS(snpcov, paste0(outdir, '/snps.rds'))
  saveRDS(gt, paste0(outdir, '/gt.rds'))
  message('dumped out files and finished')
  
}


#' @name simulate_coverage
#' @title simulate_coverage
#'
#' @description
#' simulate binned read counts given expected depth
#'
#' @details
#' Given a genomic segments with estimated tumor copy number (CN) and genomic bin size, simulate read counts per genomic bin.
#' This can be done for either a diploid genome (diploid = TRUE) or haploid genome (haploid = TRUE).
#'
#' Optionally, a multiplicative bias (representing replication timing, batch effect, etc.) can be supplied to make the coverage appear more realistic.
#' This should be supplied as a GRanges as a numeric vector in field "background" (or some other field specified in bias.field)
#' 
#' The expected coverage and read size can be tuned by setting basecov and readsize, respectively.
#'
#' By default, read counts are simulated from a Poisson distribution, but if overdispersion (numeric) is supplied will simulate from a Negative Binomial distribution.
#'
#' @param gr (GRanges, gGraph, or file containing one of these things) segment copy number, with numeric cn supplied in field "cn"
#' @param bins (GRanges) genomic bins
#' @param binsize (numeric) width of genomic bins (bp) for simulating read counts. ignored if bins if supplied. default 1e3.
#' @param diploid (logical) simulate diploid genome? default TRUE
#' @param bias (GRanges) mulitplicative bias for binned means
#' @param bias.field (character) column in bias GRanges metadata containing multiplicative bias. default "background"
#' @param basecov (numeric) expected coverage depth. default 60.
#' @param readsize (numeric) expected read size (bp). default 150.
#' @param purity (numeric) sample purity, between zero and one. default 0.95.
#' @param poisson (numeric) simulate actual integer reads? if NA will just output relative copy number.
#' @param overdispersion (numeric) if supplied, should be > 1 and will be applied to simulate from a negative binomial distribution. default NA.
#' @param normalize (logical) unit normalize resulting reads by dividing by mean? default TRUE
#
#' @return GRanges with coverage per genomic bin in field "cov" and pre-bias coverage per genomic bin in field "cov.og"
#'
#' @export
simulate_coverage = function(gr, ##needs field cn
                             bins = NULL,
                             binsize = 1000,
                             diploid = TRUE,
                             bias = NULL,
                             basecov = 60,
                             readsize = 150,
                             purity = 0.95,
                             poisson = TRUE,
                             overdispersion = NA,
                             normalize = TRUE)
{
  if (inherits(gr, 'gGraph'))
    gr = gr$nodes$gr
  #chrom.sizes = grab_chrom_sizes(genome)
  
  ## generate tiled genome
  if (is.null(bins) || !inherits(bins, "GRanges")){
    bins = gr.tile(seqlengths(gr), binsize)
  }
  
  binsize = unique(width(bins))[1]
  
  ## compute expected bin coverage
  bincov = (basecov * (readsize + binsize - 1)) / readsize
  
  ## get expected normal cn
  ncn = 1
  if (diploid) { ncn = 2 }
  
  ## compute relative CN
  GenomicRanges::mcols(bins)[, "cn"] = gUtils::gr.val(query = bins,
                                                      target = gr,
                                                      val = "cn",
                                                      mean = TRUE,
                                                      na.rm = TRUE)$cn
  
  GenomicRanges::mcols(bins)[, "cn.rel"] = (purity * GenomicRanges::mcols(bins)[, "cn"] +
                                              ncn * (1 - purity)) /
    (purity * mean(GenomicRanges::mcols(bins)[, "cn"], na.rm = TRUE) +
       ncn * (1 - purity))
  
  ## add multiplicative bias, if supplied
  #GenomicRanges::mcols(bins)[, "bias"] = 1
  if (!is.null(bias)){
    #bias = read_segs(bias, field = bias.field)
    GenomicRanges::mcols(bins)[, "bias"] = gUtils::gr.val(query = bins,
                                                                target = bias,
                                                                val = "bias",
                                                                mean = TRUE,
                                                                na.rm = TRUE)$bias
  }
  
  ## simulate from Poisson distribution
  GenomicRanges::mcols(bins)[, "cov.og"] = GenomicRanges::mcols(bins)[, "cn.rel"] * bincov #prebias
  GenomicRanges::mcols(bins)[, "cov"] = GenomicRanges::mcols(bins)[, "cn.rel"] * GenomicRanges::mcols(bins)[, "bias"] * bincov
  
  if (poisson) {
    if (is.na(overdispersion)){
      GenomicRanges::mcols(bins)[, "cov"] = stats::rpois(n = length(bins), lambda = GenomicRanges::mcols(bins)[, "cov"])
      GenomicRanges::mcols(bins)[, "cov.og"] = stats::rpois(n = length(bins), lambda = GenomicRanges::mcols(bins)[, "cov.og"])
    } else {
      GenomicRanges::mcols(bins)[, "cov"] = stats::rbinom(n = length(bins), mu =  GenomicRanges::mcols(bins)[, "cov"], size = 1 / overdispersion)
      GenomicRanges::mcols(bins)[, "cov.og"] = stats::rbinom(n = length(bins), mu =  GenomicRanges::mcols(bins)[, "cov.og"], size = 1 / overdispersion)
    }
  }
    
  ## normalize if desired
  if (normalize) {
    GenomicRanges::mcols(bins)[, "cov"] = GenomicRanges::mcols(bins)[, "cov"]  / mean(GenomicRanges::mcols(bins)[, "cov"], na.rm = TRUE)
  }
  return(bins)
}

#' @name .simcov
#'
#' @description
#'
#' @param gr input gRanges on which to simulate coverage (requires field 'cn')
#' @param binsize
#' @param bins 
#' @param bias granges of biases for given regions, first field is assumed to be the bias amount
#' @param diploid 
#' @param readsize 
#' @param basecov 
#' @param overdispersion 
#' @param purity target purity
#' @param poisson add poisson noise? default == T
#' @param normalize 
.simcov = function(gr, ## needs field 'cn'
                   binsize = 1e3,
                   bins = gr.tile(seqlengths(gr), binsize),
                   bias = NULL, 
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
  
  
  if (!is.null(bias)) {
    values(bias)[[1]] = values(bias)[[1]]/mean(values(bias)[[1]], na.rm = TRUE)
    bins$cov.og = bins$cov
    ##bins$cov = bins$cov * values(bins[, c()] %$% bias[, 1])[[1]] ## made this immune to NA for targeted seq (Addy)
    bins$cov = bins$cov * values(gr.val(bins[, c()], bias[, 1], na.rm = T, val = names(values(bias[, 1]))))[[1]]
  }
  
  return(bins)
}

#' @name .separate
#' @description 
#'
#' @param dt 
#' @param old 
#' @param new 
#' @param sep 
#' @param perl 
.separate = function(dt, old, new, sep, perl = FALSE)
{
  newdt = dt[[old]] %>% 
    as.character %>% 
    strsplit(split = sep, perl = perl) %>% 
    do.call(rbind, .) %>% 
    as.data.table %>% 
    setnames(new)
  dt[[old]] = NULL
  dt = cbind(dt, newdt)
  return(dt)
}
