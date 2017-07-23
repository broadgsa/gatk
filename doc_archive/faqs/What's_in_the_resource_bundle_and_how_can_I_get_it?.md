## What's in the resource bundle and how can I get it?

http://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it

<p><strong>NOTE: we recently made some changes to the bundle on the FTP server; see the <a href="https://software.broadinstitute.org/gatk/download/bundle">Resource Bundle</a> page for details. In a nutshell: minor directory structure changes, and Hg38 bundle now mirrors the cloud version.</strong></p>
<hr />
<h3>1. Accessing the bundle</h3>
<p>See the <a href="https://software.broadinstitute.org/gatk/download/bundle">Resource Bundle</a> page. In a nutshell, there's a Google Cloud bucket and an FTP server. The cloud bucket only has Hg38 resources; the resources for other builds are currently only available through the FTP server. Let us know if you want them on the Cloud too. </p>
<hr />
<h3>2. Grch38/Hg38 Resources: the soon-to-be Standard Set</h3>
<p>This contains all the resource files needed for Best Practices short variant discovery in whole-genome sequencing data (WGS). Exome files and itemized resource list coming soon(ish). </p>
<hr />
<h4>All resources below this are available only on the FTP server, not on the cloud.</h4>
<hr />
<h3>3. b37 Resources: the Standard Data Set pending completion of the Hg38 bundle</h3>
<ul>
<li>Reference sequence (standard 1000 Genomes fasta) along with fai and dict files</li>
<li>dbSNP in VCF.  This includes two files:
<ul>
<li>A recent dbSNP release (build 138)</li>
<li>This file subsetted to only sites discovered in or before dbSNPBuildID 129, which excludes the impact of the 1000 Genomes project and is useful for evaluation of dbSNP rate and Ti/Tv values at novel sites.</li>
</ul></li>
<li>HapMap genotypes and sites VCFs</li>
<li>OMNI 2.5 genotypes for 1000 Genomes samples, as well as sites, VCF </li>
<li>The current best set of known indels to be used for local realignment (note that we don't use dbSNP for this anymore); use both files:
<ul>
<li>1000G_phase1.indels.b37.vcf (currently from the 1000 Genomes Phase I indel calls)</li>
<li>Mills_and_1000G_gold_standard.indels.b37.sites.vcf</li>
</ul></li>
<li>The latest set from 1000G phase 3 (v4) for genotype refinement: 1000G_phase3_v4_20130502.sites.vcf</li>
<li>A large-scale standard single sample BAM file for testing:
<ul>
<li>NA12878.HiSeq.WGS.bwa.cleaned.recal.b37.20.bam containing ~64x reads of NA12878 on chromosome 20</li>
<li>A callset produced by running UnifiedGenotyper on the dataset above. Note that this resource is out of date and does not represent the results of our Best Practices. This will be updated in the near future.</li>
</ul></li>
<li>The Broad's custom exome targets list: Broad.human.exome.b37.interval_list (note that you should always use the exome targets list that is appropriate for your data, which typically depends on the prep kit that was used, and should be available from the kit manufacturer's website)</li>
</ul>
<p>Additionally, these files all have supplementary indices, statistics, and other QC data available.</p>
<hr />
<h3>4. hg19 Resources: lifted over from b37</h3>
<p>Includes the UCSC-style hg19 reference along with all lifted over VCF files.</p>
<hr />
<h3>5. hg18 Resources: lifted over from b37</h3>
<p>Includes the UCSC-style hg18 reference along with all lifted over VCF files. The refGene track and BAM files are not available. We only provide data files for this genome-build that can be lifted over &quot;easily&quot; from our master b37 repository.  Sorry for whatever inconvenience that this might cause.</p>
<p>Also includes a chain file to lift over to b37.</p>
<hr />
<h3>6. b36 Resources: lifted over from b37</h3>
<p>Includes the 1000 Genomes pilot b36 formatted reference sequence (human_b36_both.fasta) along with all lifted over VCF files. The refGene track and BAM files are not available.  We only provide data files for this genome-build that can be lifted over &quot;easily&quot; from our master b37 repository.  Sorry for whatever inconvenience that this might cause.</p>
<p>Also includes a chain file to lift over to b37.</p>