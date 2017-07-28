## How should I pre-process data from multiplexed sequencing and multi-library designs?

http://gatkforums.broadinstitute.org/gatk/discussion/3060/how-should-i-pre-process-data-from-multiplexed-sequencing-and-multi-library-designs

<p>Our Best Practices pre-processing documentation assumes a simple experimental design in which you have one set of input sequence files (forward/reverse or interleaved FASTQ, or unmapped uBAM) per sample, and you run each step of the pre-processing workflow separately for each sample, resulting in one BAM file per sample at the end of this phase. </p>
<p>However, if you are generating multiple libraries for each sample, and/or multiplexing samples within and/or across sequencing lanes, the data must be de-multiplexed before pre-processing, typically resulting in multiple sets of FASTQ files per sample all of which should have distinct <a href="https://www.broadinstitute.org/gatk/guide/article?id=6472">read group</a> IDs (RGID). </p>
<p>At that point there are several different valid strategies for implementing the pre-processing workflow. Here at the Broad Institute, we run the initial steps of the pre-processing workflow (mapping, sorting and marking duplicates) separately on each individual read group. Then we merge the data to produce a single BAM file for each sample (aggregation); this is done by re-running Mark Duplicates, this time on all read group BAM files for a sample at the same time. Then we run Indel Realignment and Base Recalibration on the aggregated per-sample BAM files. See the worked-out example below and <a href="https://www.broadinstitute.org/gatk/events/slides/1506/GATKwr8-A-3-GATK_Best_Practices_and_Broad_pipelines.pdf">this presentation</a> for more details.</p>
<p><em>Note that there are many possible ways to achieve a similar result; here we present the way we think gives the best combination of efficiency and quality. This assumes that you are dealing with one or more samples, and each of them was sequenced on one or more lanes.</em></p>
<h3>Example</h3>
<p>Let's say we have this example data (assuming interleaved FASTQs containing both forward and reverse reads) for two sample libraries, <em>sampleA</em> and <em>sampleB</em>, which were each sequenced on two lanes, <em>lane1</em> and <em>lane2</em>:</p>
<ul>
<li>sampleA_lane1.fq</li>
<li>sampleA_lane2.fq</li>
<li>sampleB_lane1.fq</li>
<li>sampleB_lane2.fq</li>
</ul>
<p>These will each be identified as separate read groups A1, A2, B1 and B2. If we had multiple libraries per sample, we would further distinguish them (eg sampleA_lib1_lane1.fq leading to read group A11, sampleA_lib2_lane1.fq leading to read group A21 and so on).</p>
<h4>1. Run initial steps per-readgroup once</h4>
<p>Assuming that you received one FASTQ file per sample library, per lane of sequence data (which amounts to a <a href="https://www.broadinstitute.org/gatk/guide/article?id=6472">read group</a>), run each file through mapping and  sorting. During the mapping step you assign read group information, which will be very important in the next steps so be sure to do it correctly. See the <a href="https://www.broadinstitute.org/gatk/guide/article?id=6472">read groups</a> dictionary entry for guidance. </p>
<p>The example data becomes:</p>
<ul>
<li>sampleA_rgA1.bam</li>
<li>sampleA_rgA2.bam</li>
<li>sampleB_rgB1.bam</li>
<li>sampleB_rgB2.bam</li>
</ul>
<p>At this point we mark duplicates in each read group BAM file (dedup), which allows us to estimate the complexity of the corresponding library of origin as a quality control step. This step is optional. </p>
<p>The example data becomes:</p>
<ul>
<li>sampleA_rgA1.dedup.bam</li>
<li>sampleA_rgA2.dedup.bam</li>
<li>sampleB_rgB1.dedup.bam</li>
<li>sampleB_rgB2.dedup.bam</li>
</ul>
<p>Technically this first run of marking duplicates is not necessary because we will run it again per-sample, and that per-sample marking would be enough to achieve the desired result. To reiterate, we only do this round of marking duplicates for QC purposes. </p>
<h4>2. Merge read groups and mark duplicates per sample (aggregation + dedup)</h4>
<p>Once you have pre-processed each read group individually, you merge read groups belonging to the same sample into a single BAM file. You can do this as a standalone step, bur for the sake of efficiency we combine this with the per-readgroup duplicate marking step (it's simply a matter of passing the multiple inputs to MarkDuplicates in a single command). </p>
<p>The example data becomes:</p>
<ul>
<li>sampleA.merged.dedup.bam</li>
<li>sampleB.merged.dedup.bam</li>
</ul>
<p>To be clear, this is the round of marking duplicates that matters. It eliminates PCR duplicates (arising from library preparation) across all lanes in addition to optical duplicates (which are by definition only per-lane). </p>
<h4>3. Remaining per-sample pre-processing</h4>
<p>Then you run indel realignment (optional) and base recalibration (BQSR). </p>
<p>The example data becomes:</p>
<ul>
<li>sample1.merged.dedup.(realn).recal.bam</li>
<li>sample2.merged.dedup.(realn).recal.bam</li>
</ul>
<p>Realigning around indels per-sample leads to consistent alignments across all lanes within a sample. This step is only necessary if you will be using a locus-based variant caller like MuTect 1 or UnifiedGenotyper (for legacy reasons). If you will be using HaplotypeCaller or MuTect2, you do not need to perform indel realignment. </p>
<p>Base recalibration will be applied per-read group if you assigned appropriate read group information in your data. BaseRecalibrator distinguishes read groups by RGID, or RGPU if it is available (PU takes precedence over ID). This will identify separate read groups (distinguishing both lanes and libraries) as such even if they are in the same BAM file, and it will always process them separately -- as long as the read groups are identified correctly of course. There would be no sense in trying to recalibrate across lanes, since the purpose of this processing step is to compensate for the errors made by the machine during sequencing, and the lane is the base unit of the sequencing machine (assuming the equipment is Illumina HiSeq or similar technology). </p>
<p><em>People often ask also if it's worth the trouble to try realigning across all samples in a cohort. The answer is almost always no, unless you have very shallow coverage. The problem is that while it would be lovely to ensure consistent alignments around indels across all samples, the computational cost gets too ridiculous too fast. That being said, for contrastive calling projects -- such as cancer tumor/normals -- we do recommend realigning both the tumor and the normal together in general to avoid slight alignment differences between the two tissue types.</em></p>