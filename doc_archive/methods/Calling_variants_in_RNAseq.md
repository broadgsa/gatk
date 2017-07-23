## Calling variants in RNAseq

http://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq

<h3>Overview</h3>
<p>This document describes the details of the GATK Best Practices workflow for SNP and indel calling on RNAseq data.</p>
<p><strong>Please note that any command lines are only given as example of how the tools can be run. You should always make sure you understand what is being done at each step and whether the values are appropriate for your data. To that effect, you can find more guidance <a href="http://www.broadinstitute.org/gatk/guide/best-practices?bpm=RNAseq">here</a>.</strong> </p>
<p><a href='https://us.v-cdn.net/5019796/uploads/FileUpload/fa/e60ecf89bd1b2645d9fce68ccf3919.png'><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/fa/e60ecf89bd1b2645d9fce68ccf3919.png" /></a></p>
<p>In brief, the key modifications made to the DNAseq Best Practices focus on handling splice junctions correctly, which involves specific mapping and pre-processing procedures, as well as some new functionality in the HaplotypeCaller. Here is a detailed overview:</p>
<p><a href='https://us.v-cdn.net/5019796/uploads/FileUpload/c9/ac46784be39f31fa976b5ac944de17.png'><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/c9/ac46784be39f31fa976b5ac944de17.png" /></a></p>
<h3>Caveats</h3>
<p>Please keep in mind that our DNA-focused Best Practices were developed over several years of thorough experimentation, and are continuously updated as new observations come to light and the analysis methods improve. We have been working with RNAseq for a somewhat shorter time, so there are many aspects that we still need to examine in more detail before we can be fully confident that we are doing the best possible thing. </p>
<p>We know that the current recommended pipeline is producing both false positives (wrong variant calls) and false negatives (missed variants) errors. While some of those errors are inevitable in any pipeline, others are errors that we can and will address in future versions of the pipeline. A few examples of such errors are given in this article as well as our ideas for fixing them in the future.  </p>
<p>We will be improving these recommendations progressively as we go, and we hope that the research community will help us by providing feedback of their experiences applying our recommendations to their data.</p>
<hr />
<h3>The workflow</h3>
<h4>1. Mapping to the reference</h4>
<p>The first major difference relative to the DNAseq Best Practices is the mapping step. For DNA-seq, we recommend BWA. For RNA-seq, we evaluated all the major software packages that are specialized in RNAseq alignment, and we found that we were able to achieve the highest sensitivity to both SNPs and, importantly, indels, using STAR aligner. Specifically, we use the STAR 2-pass method which was described in a recent publication (see page 43 of the Supplemental text of the Pär G Engström et al. paper referenced below for full protocol details -- we used the suggested protocol with the default parameters). In brief, in the STAR 2-pass approach, splice junctions detected in a first alignment run are used to guide the final alignment.</p>
<p>Here is a walkthrough of the STAR 2-pass alignment steps:</p>
<p>1) STAR uses genome index files that must be saved in unique directories. The human genome index was built from the FASTA file hg19.fa as follows:</p>
<pre><code class="pre_md">genomeDir=/path/to/hg19
mkdir $genomeDir
STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles hg19.fa\  --runThreadN &lt;n&gt;</code class="pre_md"></pre>
<p>2) Alignment jobs were executed as follows:</p>
<pre><code class="pre_md">runDir=/path/to/1pass
mkdir $runDir
cd $runDir
STAR --genomeDir $genomeDir --readFilesIn mate1.fq mate2.fq --runThreadN &lt;n&gt;</code class="pre_md"></pre>
<p>3) For the 2-pass STAR, a new index is then created using splice junction information contained in the file SJ.out.tab from the first pass:</p>
<pre><code class="pre_md">genomeDir=/path/to/hg19_2pass
mkdir $genomeDir
STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles hg19.fa \
    --sjdbFileChrStartEnd /path/to/1pass/SJ.out.tab --sjdbOverhang 75 --runThreadN &lt;n&gt;</code class="pre_md"></pre>
<p>4) The resulting index is then used to produce the final alignments as follows:</p>
<pre><code class="pre_md">runDir=/path/to/2pass
mkdir $runDir
cd $runDir
STAR --genomeDir $genomeDir --readFilesIn mate1.fq mate2.fq --runThreadN &lt;n&gt;</code class="pre_md"></pre>
<h4>2. Add read groups, sort, mark duplicates, and create index</h4>
<p>The above step produces a SAM file, which we then put through the usual Picard processing steps: adding read group information, sorting, marking duplicates and indexing. </p>
<pre><code class="pre_md">java -jar picard.jar AddOrReplaceReadGroups I=star_output.sam O=rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample 

java -jar picard.jar MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics </code class="pre_md"></pre>
<h4>3. Split'N'Trim and reassign mapping qualities</h4>
<p>Next, we use a new GATK tool called SplitNCigarReads developed specially for RNAseq, which splits reads into exon segments (getting rid of Ns but maintaining grouping information) and hard-clip any sequences overhanging into the intronic regions. </p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/60/d400461fa1673ce50603487714be76.png" />
<p>In the future we plan to integrate this into the GATK engine so that it will be done automatically where appropriate, but for now it needs to be run as a separate step.</p>
<p>At this step we also add one important tweak: we need to reassign mapping qualities, because STAR assigns good alignments a MAPQ of 255 (which technically means “unknown” and is therefore meaningless to GATK). So we use the GATK’s <a href="http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_filters_ReassignOneMappingQualityFilter.html">ReassignOneMappingQuality</a> read filter to reassign all good alignments to the default value of 60. This is not ideal, and we hope that in the future RNAseq mappers will emit meaningful quality scores, but in the meantime this is the best we can do. In practice we do this by adding the ReassignOneMappingQuality read filter to the splitter command.</p>
<p>Finally, be sure to specify that reads with N cigars should be allowed. This is currently still classified as an &quot;unsafe&quot; option, but this classification will change to reflect the fact that this is now a supported option for RNAseq processing.</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R ref.fasta -I dedupped.bam -o split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS</code class="pre_md"></pre>
<h4>4. Indel Realignment (optional)</h4>
<p>After the splitting step, we resume our regularly scheduled programming... to some extent. We have found that performing realignment around indels can help rescue a few indels that would otherwise be missed, but to be honest the effect is marginal. So while it can’t hurt to do it, we only recommend performing the realignment step if you have compute and time to spare (or if it’s important not to miss any potential indels).</p>
<h4>5. Base Recalibration</h4>
<p>We do recommend running base recalibration (BQSR). Even though the effect is also marginal when applied to good quality data, it can absolutely save your butt in cases where the qualities have systematic error modes. </p>
<p>Both steps 4 and 5  are run as described for DNAseq (with the same known sites resource files), without any special arguments. Finally, please note that you should NOT run ReduceReads on your RNAseq data. The ReduceReads tool will no longer be available in GATK 3.0.</p>
<h4>6. Variant calling</h4>
<p>Finally, we have arrived at the variant calling step! Here, we recommend using HaplotypeCaller because it is performing much better in our hands than UnifiedGenotyper (our tests show that UG was able to call less than 50% of the true positive indels that HC calls). We have added some functionality to the variant calling code which will intelligently take into account the information about intron-exon split regions that is embedded in the BAM file by SplitNCigarReads. In brief, the new code will perform “dangling head merging” operations and avoid using soft-clipped bases (this is a temporary solution) as necessary to minimize false positive and false negative calls. To invoke this new functionality, just add <code>-dontUseSoftClippedBases</code> to your regular HC command line. Note that the <code>-recoverDanglingHeads</code> argument which was previously required is no longer necessary as that behavior is now enabled by default in HaplotypeCaller. Also, we found that we get better results if we set the minimum phred-scaled confidence threshold for calling variants 20, but you can lower this to increase sensitivity if needed.</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R ref.fasta -I input.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o output.vcf</code class="pre_md"></pre>
<h4>7. Variant filtering</h4>
<p>To filter the resulting callset, you will need to apply hard filters, as we do not yet have the RNAseq training/truth resources that would be needed to run variant recalibration (VQSR). </p>
<p>We recommend that you filter clusters of at least 3 SNPs that are within a window of 35 bases between them by adding <code>-window 35 -cluster 3</code> to your command. This filter recommendation is specific for RNA-seq data. </p>
<p>As in DNA-seq, we recommend filtering based on Fisher Strand values (FS &gt; 30.0) and Qual By Depth values (QD &lt; 2.0).</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T VariantFiltration -R hg_19.fasta -V input.vcf -window 35 -cluster 3 -filterName FS -filter "FS &gt; 30.0" -filterName QD -filter "QD &lt; 2.0" -o output.vcf </code class="pre_md"></pre>
<p>Please note that we selected these hard filtering values in attempting to optimize both high sensitivity and specificity together. By applying the hard filters, some real sites will get filtered. This is a tradeoff that each analyst should consider based on his/her own project.  If you care more about sensitivity and are willing to tolerate more false positives calls, you can choose not to filter at all (or to use less restrictive thresholds).     </p>
<p>An example of filtered (SNPs cluster filter) and unfiltered false variant calls: </p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/08/a745a2deb39ff5d9ab238f0cd9c58b.png" />
<p>An example of true variants that were filtered (false negatives). As explained in text, there is a tradeoff that comes with applying filters:</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/09/5c95bdbe402f2c903420ad4d599bba.png" />
<hr />
<h3>Known issues</h3>
<p>There are a few known issues; one is that the allelic ratio is problematic. In many heterozygous sites, even if we can see in the RNAseq data both alleles that are present in the DNA, the ratio between the number of reads with the different alleles is far from 0.5, and thus the HaplotypeCaller (or any caller that expects a diploid genome) will miss that call. A DNA-aware mode of the caller might be able to fix such cases (which may be candidates also for downstream analysis of allele specific expression). </p>
<p>Although our new tool (splitNCigarReads) cleans many false positive calls that are caused by splicing inaccuracies by the aligners, we still call some false variants for that same reason, as can be seen in the example below. Some of those errors might be fixed in future versions of the pipeline with more sophisticated filters, with another realignment step in those regions, or by making the caller aware of splice positions.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/98/3ac331af5f6b24ead9acd56a546e22.png" />
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/30/35d57e11ae6f638f17043c5d2f09ef.png" />
<p>As stated previously, we will continue to improve the tools and process over time. We have plans to improve the splitting/clipping functionalities, improve true positive and minimize false positive rates, as well as developing statistical filtering (i.e. variant recalibration) recommendations. </p>
<p>We also plan to add functionality to process DNAseq and RNAseq data from the same samples simultaneously, in order to facilitate analyses of post-transcriptional processes. Future extensions to the HaplotypeCaller will provide this functionality, which will require both DNAseq and RNAseq in order to produce the best results. Finally, we are also looking at solutions for measuring differential expression of alleles. </p>
<hr />
<p><em>[1] Pär G Engström et al. “Systematic evaluation of spliced alignment programs for RNA-seq data”. Nature Methods, 2013</em></p>
<hr />
<h4>NOTE: Questions about this document that were posted before June 2014 have been moved to this archival thread: <a href="http://gatkforums.broadinstitute.org/discussion/4709/questions-about-the-rnaseq-variant-discovery-workflow">http://gatkforums.broadinstitute.org/discussion/4709/questions-about-the-rnaseq-variant-discovery-workflow</a></h4>