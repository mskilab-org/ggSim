~~~

   ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄       ▄▄  ▄ 
  ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░▌     ▐░░▌▐░▌
  ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀  ▀▀▀▀█░█▀▀▀▀ ▐░▌░▌   ▐░▐░▌▐░▌
  ▐░▌          ▐░▌          ▐░▌               ▐░▌     ▐░▌▐░▌ ▐░▌▐░▌▐░▌
  ▐░▌ ▄▄▄▄▄▄▄▄ ▐░▌ ▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄▄▄      ▐░▌     ▐░▌ ▐░▐░▌ ▐░▌▐░▌
  ▐░▌▐░░░░░░░░▌▐░▌▐░░░░░░░░▌▐░░░░░░░░░░░▌     ▐░▌     ▐░▌  ▐░▌  ▐░▌▐░▌
  ▐░▌ ▀▀▀▀▀▀█░▌▐░▌ ▀▀▀▀▀▀█░▌ ▀▀▀▀▀▀▀▀▀█░▌     ▐░▌     ▐░▌   ▀   ▐░▌▐░▌
  ▐░▌       ▐░▌▐░▌       ▐░▌          ▐░▌     ▐░▌     ▐░▌       ▐░▌ ▀ 
  ▐░█▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄█░▌ ▄▄▄▄▄▄▄▄▄█░▌ ▄▄▄▄█░█▄▄▄▄ ▐░▌       ▐░▌ ▄ 
  ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░▌       ▐░▌▐░▌
   ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀         ▀  ▀ 
                                                                    
~~~

---

## Robust simulation of whole-genome gGraphs featuring
 * Coverage of tumor and normal tracks
 * SNP phasing
 * Junction phasing
 * Coverage of SNPs

## <font color=black> Tutorial </font>

`ggsim` is the main function in this package to simulate a genome. Given a set of junctions and short-nucleotide polymorphisms, `ggsim` will create robust sex-informed, junction-balanced phased and unphased gGraphs with corresponding coverages that match the user's input purity and ploidy.   

The essential parameters to this function are `junctions`, `vcf`, `bias`, and `nbias`, and must be supplied. `junctions` is a Junctions object, `vcf` defines the SNP profile for the sim genome, and `bias`/`nbias` represent the normal coverage vectors multiplied to the tumor/normal coverage, respectively, to simulate real-world fluctuations in read depth. More details about function parameters are outlined in the table below.

<table style="border: 1px solid black; border-collapse: collapse;">
  <tbody>
    <tr>
      <th style="border: 1px solid black; padding: 5px;">Parameter</th>
      <th style="border: 1px solid black; padding: 5px;">Default value</th>
      <th style="border: 1px solid black; padding: 5px;">Description/notes</th>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>junctions</code></td>
      <td style="border: 1px solid black; padding: 5px;"></td>
      <td style="border: 1px solid black; padding: 5px;">Junctions to add to gGraph as a GRangesList</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>vcf</code></td>
      <td style="border: 1px solid black; padding: 5px;"></td>
      <td style="border: 1px solid black; padding: 5px;">Phased VCF of germline heterozygous SNPs. Can use any pileup of a genome or a [Platinum Genome from Illumina](https://github.com/Illumina/PlatinumGenomes).<br><br> <b>NOTE:</b> the sex of this input determines the sex of the simulated genome. Presence/absence of heterozygous SNPs will define genome as F/M, with subsequent effects on the defined CN/haplotyping of sex chromosomes. </td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>bias</code></td>
      <td style="border: 1px solid black; padding: 5px;"></td>
      <td style="border: 1px solid black; padding: 5px;">.rds of binned read depth bias for tumor sample e.g. read depth for a random normal sample</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>nbias</code></td>
      <td style="border: 1px solid black; padding: 5px;"></td>
      <td style="border: 1px solid black; padding: 5px;">.rds of binned read depth bias for normal sample e.g. read depth for a random normal sample</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>snps</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>NULL</code></td>
      <td style="border: 1px solid black; padding: 5px;">Optional comprehensive VCF of reference snps e.g. hapmap</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>unmappable</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>NULL</code></td>
      <td style="border: 1px solid black; padding: 5px;">Optional .rds of GRanges of CN unmappable regions</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>coverage</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>60</code></td>
      <td style="border: 1px solid black; padding: 5px;">Target tumor base coverage</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>ncoverage</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>40</code></td>
      <td style="border: 1px solid black; padding: 5px;">Target normal base coverage</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>alpha</code>\
      (purity)</td>
      <td style="border: 1px solid black; padding: 5px;"><code>1.00</code></td>
      <td style="border: 1px solid black; padding: 5px;">Target purity</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>tau</code>\
      (ploidy)</td>
      <td style="border: 1px solid black; padding: 5px;"><code>1.00</code></td>
      <td style="border: 1px solid black; padding: 5px;">Target ploidy</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>poisson</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>TRUE</code></td>
      <td style="border: 1px solid black; padding: 5px;">Add shot noise to read depth?</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>numbreaks</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>10</code></td>
      <td style="border: 1px solid black; padding: 5px;">Number of additional breaks to add in CN-unmappable regions</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>width</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>1000</code></td>
      <td style="border: 1px solid black; padding: 5px;">Bin width of read depth</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>cnloh</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>FALSE</code></td>
      <td style="border: 1px solid black; padding: 5px;">Add a copy-neutral loss of heterozygosity edge?</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>standard.chr</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>c(1:22, "X", "Y")</code></td>
      <td style="border: 1px solid black; padding: 5px;">Defaults to human chromosomes</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>outdir</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>./</code></td>
      <td style="border: 1px solid black; padding: 5px;">Path to save</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;"><code>par.path</code></td>
      <td style="border: 1px solid black; padding: 5px;"><code>system.file("extdata", "PAR_hg19.rds", package = 'ggSim')</code></td>
      <td style="border: 1px solid black; padding: 5px;">GRanges identifying pseudoautosomal regions in X, Y chromosomes. This is used to agnosticize the sex of the <code>bias</code>/<code>nbias</code> vectors.</td>
    </tr>
  </tbody>
</table>


