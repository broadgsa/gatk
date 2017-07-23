## Errors about BAM or VCF files not being ordered properly

http://gatkforums.broadinstitute.org/gatk/discussion/58/errors-about-bam-or-vcf-files-not-being-ordered-properly

<h3>This article has been deprecated</h3>
<h4>For a more recent version please see <a href="https://www.broadinstitute.org/gatk/guide/article?id=1328">https://www.broadinstitute.org/gatk/guide/article?id=1328</a></h4>
<hr />
<p>This error occurs when for example, a collaborator gives you a BAM that's derived from what was originally the same reference as you are using, but for whatever reason the contigs are not sorted in the same order .The GATK can be particular about the <a href="http://www.broadinstitute.org/gatk/guide/article?id=1204">ordering of a BAM file</a> so it will fail with an error in this case. </p>
<p>So what do you do? You use a Picard tool called ReorderSam to, well, reorder your BAM file. </p>
<p>Here's an example usage where we reorder a BAM file that was sorted lexicographically so that the output will be another BAM, but this time sorted karyotypically : </p>
<pre><code class="pre_md">java -jar picard.jar ReorderSam \
    I= lexicographic.bam \
    O= kayrotypic.bam \
    REFERENCE= Homo_sapiens_assembly18.kayrotypic.fasta</code class="pre_md"></pre>
<p>This tool requires you have a correctly sorted version of the reference sequence you used to align your reads.  Be aware that this tool will drop reads that don't have equivalent contigs in the new reference (potentially bad, but maybe not). If contigs have the same name in the bam and the new reference, this tool assumes that the alignment of the read in the new BAM is the same.  This is not a liftover tool!</p>
<p>This tool is part of the <a href="https://broadinstitute.github.io/picard/command-line-overview.html#ReorderSam">Picard package</a>.</p>