## Merging batched call sets - RETIRED

http://gatkforums.broadinstitute.org/gatk/discussion/46/merging-batched-call-sets-retired

<h3>This procedure is deprecated since it is no longer necessary and goes against our Best Practices recommendations. For calling variants on multiple samples, use the <a href="https://www.broadinstitute.org/gatk/guide/best-practices">Best Practices workflow</a> for performing variant discovery using HaplotypeCaller.</h3>
<hr />
<h3>Introduction</h3>
<p>Three-stage procedure:</p>
<ul>
<li>
<p>Create a master set of sites from your N batch VCFs that you want to genotype in all samples.  At this stage you need to determine how you want to resolve disagreements among the VCFs.  This is your master sites VCF.</p>
</li>
<li>
<p>Take the master sites VCF and genotype each sample BAM file at these sites</p>
</li>
<li>(Optionally) Merge the single sample VCFs into a master VCF file</li>
</ul>
<h3>Creating the master set of sites: SNPs and Indels</h3>
<p>The first step of batch merging is to create a master set of sites that you want to genotype in all samples.  To make this problem concrete, suppose I have two VCF files:</p>
<p>Batch 1:</p>
<pre><code class="pre_md">##fileformat=VCFv4.0
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12891 
20      9999996     .       A       ATC     .       PASS    .       GT:GQ   0/1:30
20      10000000        .       T       G       .       PASS    .       GT:GQ   0/1:30
20      10000117        .       C       T       .       FAIL    .       GT:GQ   0/1:30
20      10000211        .       C       T       .       PASS    .       GT:GQ   0/1:30
20      10001436        .       A       AGG     .       PASS    .       GT:GQ   1/1:30</code class="pre_md"></pre>
<p>Batch 2:</p>
<pre><code class="pre_md">##fileformat=VCFv4.0
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12878
20      9999996     .       A       ATC     .       PASS    .       GT:GQ   0/1:30
20      10000117        .       C       T       .       FAIL    .       GT:GQ   0/1:30
20      10000211        .       C       T       .       FAIL    .       GT:GQ   0/1:30
20      10000598        .       T       A       .       PASS    .       GT:GQ   1/1:30
20      10001436        .       A       AGGCT   .       PASS    .       GT:GQ   1/1:30</code class="pre_md"></pre>
<p>In order to merge these batches, I need to make a variety of bookkeeping and filtering decisions, as outlined in the merged VCF below: </p>
<p>Master VCF:</p>
<pre><code class="pre_md">20      9999996     .       A       ATC     .       PASS    .       GT:GQ   0/1:30  [pass in both]
20      10000000        .       T       G       .       PASS    .       GT:GQ   0/1:30  [only in batch 1]
20      10000117        .       C       T       .       FAIL    .       GT:GQ   0/1:30  [fail in both]
20      10000211        .       C       T       .       FAIL    .       GT:GQ   0/1:30  [pass in 1, fail in 2, choice in unclear]
20      10000598        .       T       A       .       PASS    .       GT:GQ   1/1:30  [only in batch 2]
20      10001436        .       A       AGGCT   .       PASS    .       GT:GQ   1/1:30  [A/AGG in batch 1, A/AGGCT in batch 2, including this site may be problematic]</code class="pre_md"></pre>
<p>These issues fall into the following categories:</p>
<ul>
<li>For sites present in all VCFs (20:9999996 above), the alleles agree, and each site PASS is pass, this site can obviously be considered &quot;PASS&quot; in the master VCF</li>
<li>Some sites may be PASS in one batch, but absent in others (20:10000000 and 20:10000598), which occurs when the site is polymorphic in one batch but all samples are reference or no-called in the other batch</li>
<li>Similarly, sites that are fail in all batches in which they occur can be safely filtered out, or included as failing filters in the master VCF (20:10000117)</li>
</ul>
<p>There are two difficult situations that must be addressed by the needs of the project merging batches:</p>
<ul>
<li>Some sites may be PASS in some batches but FAIL in others.  This might indicate that either:</li>
<li>The site is actually truly polymorphic, but due to limited coverage, poor sequencing, or other issues it is flag as unreliable in some batches.  In these cases, it makes sense to include the site</li>
<li>The site is actually a common machine artifact, but just happened to escape standard filtering in a few batches.  In these cases, you would obviously like to filter out the site</li>
<li>Even more complicated, it is possible that in the PASS batches you have found a reliable allele (C/T, for example) while in others there's no alt allele but actually a low-frequency error, which is flagged as failing.  Ideally, here you could filter out the failing allele from the FAIL batches, and keep the pass ones</li>
<li>Some sites may have multiple segregating alleles in each batch.  Such sites are often errors, but in some cases may be actual multi-allelic sites, in particular for indels.</li>
</ul>
<p>Unfortunately, we cannot determine which is actually the correct choice, especially given the goals of the project.  We leave it up the project bioinformatician to handle these cases when creating the master VCF.  We are hopeful that at some point in the future we'll have a consensus approach to handle such merging, but until then this will be a manual process.</p>
<p>The GATK tool <a href="http://www.broadinstitute.org/gatk/guide/article?id=53">CombineVariants</a> can be used to merge multiple VCF files, and parameter choices will allow you to handle some of the above issues.  With tools like <a href="http://www.broadinstitute.org/gatk/guide/article?id=54">SelectVariants</a>  one can slice-and-dice the merged VCFs to handle these complexities as appropriate for your project's needs.  For example, the above master merge can be produced with the following CombineVariants:</p>
<pre><code class="pre_md">java -jar dist/GenomeAnalysisTK.jar \
-T CombineVariants \
-R human_g1k_v37.fasta \
-V:one,VCF combine.1.vcf -V:two,VCF combine.2.vcf \
--sites_only \
-minimalVCF \
-o master.vcf</code class="pre_md"></pre>
<p>producing the following VCF:</p>
<pre><code class="pre_md">##fileformat=VCFv4.0
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
20      9999996     .       A       ACT         .       PASS    set=Intersection
20      10000000        .       T       G           .   PASS    set=one
20      10000117        .       C       T           .       FAIL    set=FilteredInAll
20      10000211        .       C       T           .       PASS    set=filterIntwo-one
20      10000598        .       T       A           .       PASS    set=two
20      10001436        .       A       AGG,AGGCT       .       PASS    set=Intersection</code class="pre_md"></pre>
<h3>Genotyping your samples at these sites</h3>
<p>Having created the master set of sites to genotype, along with their alleles, as in the previous section, you now use the <a href="http://www.broadinstitute.org/gatk/guide/article?id=1237">UnifiedGenotyper</a> to genotype each sample independently at the master set of sites.  This GENOTYPE_GIVEN_ALLELES mode of the UnifiedGenotyper will jump into the sample BAM file, and calculate the genotype and genotype likelihoods of the sample at the site for each of the genotypes available for the REF and ALT alleles.  For example, for site 10000211, the UnifiedGenotyper would evaluate the likelihoods of the CC, CT, and TT genotypes for the sample at this site, choose the most likely configuration, and generate a VCF record containing the genotype call and the likelihoods for the three genotype configurations.  </p>
<p>As a concrete example command line, you can genotype the master.vcf file using in the bundle sample NA12878 with the following command:</p>
<pre><code class="pre_md">java -Xmx2g -jar dist/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R bundle/b37/human_g1k_v37.fasta \
-I bundle/b37/NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.bam \
-alleles master.vcf \
-L master.vcf \
-gt_mode GENOTYPE_GIVEN_ALLELES \
-out_mode EMIT_ALL_SITES \
-stand_call_conf 0.0 \
-glm BOTH \
-G none \</code class="pre_md"></pre>
<p>The <code>-L master.vcf</code> argument tells the UG to only genotype the sites in the master file. If you don't specify this, the UG will genotype the master sites in GGA mode, but it will also genotype all other sites in the genome in regular mode. </p>
<p><code>The last item,</code>-G ` prevents the UG from computing annotations you don't need.  This command produces something like the following output:</p>
<pre><code class="pre_md">##fileformat=VCFv4.0
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12878
20      9999996     .       A       ACT         4576.19 .       .   GT:DP:GQ:PL     1/1:76:99:4576,229,0
20      10000000        .       T       G           0       .       .       GT:DP:GQ:PL     0/0:79:99:0,238,3093
20      10000211        .       C       T       857.79  .       .   GT:AD:DP:GQ:PL  0/1:28,27:55:99:888,0,870
20      10000598        .       T       A           1800.57 .       .   GT:AD:DP:GQ:PL  1/1:0,48:48:99:1834,144,0
20      10001436        .       A       AGG,AGGCT       1921.12 .       .   GT:DP:GQ:PL     0/2:49:84.06:1960,2065,0,2695,222,84</code class="pre_md"></pre>
<p>Several things should be noted here:</p>
<ul>
<li>The genotype likelihoods calculation evolves, especially for indels, the exact results of this command will change.  </li>
<li>The command will emit sites that are hom-ref in the sample at the site, but the -stand_call_conf 0.0 argument should be provided so that they aren't tagged as &quot;LowQual&quot; by the UnifiedGenotyper.</li>
<li>The filtered site 10000117 in the master.vcf is not genotyped by the UG, as it doesn't pass filters and so is considered bad by the GATK UG.  If you want to determine the genotypes for all sites, independent on filtering, you must unfilter all of your records in master.vcf, and if desired, restore the filter string for these records later.</li>
</ul>
<p>This genotyping command can be performed independently per sample, and so can be parallelized easily on a farm with one job per sample, as in the following:</p>
<pre><code class="pre_md">foreach sample in samples:
  run UnifiedGenotyper command above with -I $sample.bam -o $sample.vcf
end</code class="pre_md"></pre>
<h3>(Optional) Merging the sample VCFs together</h3>
<p>You can use a similar command for <a href="http://www.broadinstitute.org/gatk/guide/article?id=53">CombineVariants</a> above to merge back together all of your single sample genotyping runs.  Suppose all of my UnifiedGenotyper jobs have completed, and I have VCF files named sample1.vcf, sample2.vcf, to sampleN.vcf.  The single command:</p>
<pre><code class="pre_md">java -jar dist/GenomeAnalysisTK.jar -T CombineVariants -R human_g1k_v37.fasta -V:sample1 sample1.vcf -V:sample2 sample2.vcf [repeat until] -V:sampleN sampleN.vcf -o combined.vcf</code class="pre_md"></pre>
<h3>General notes</h3>
<ul>
<li>Because the GATK uses dynamic downsampling of reads, it is possible for truly marginal calls to change likelihoods from discovery (processing the BAM incrementally) vs. genotyping (jumping into the BAM).  Consequently, do not be surprised to see minor differences in the genotypes for samples from discovery and genotyping.</li>
<li>More advanced users may want to consider group several samples together for genotyping.  For example, 100 samples could be genotyped in 10 groups of 10 samples, resulting in only 10 VCF files.  Merging the 10 VCF files may be faster (or just easier to manage) than 1000 individual VCFs.</li>
<li>Sometimes, using this method, a monomorphic site within a batch will be identified as polymorphic in one or more samples within that same batch. This is because the UnifiedGenotyper applies a frequency prior to determine whether a site is likely to be monomorphic. If the site is monomorphic, it is either not output, or if EMIT_ALL_SITES is thrown, reference genotypes are output. If the site is determined to be polymorphic, genotypes are assigned greedily (as of GATK-v1.4). Calling single-sample reduces the effect of the prior, so sites which were considered monomorphic within a batch could be considered polymorphic within a sub-batch.</li>
</ul>