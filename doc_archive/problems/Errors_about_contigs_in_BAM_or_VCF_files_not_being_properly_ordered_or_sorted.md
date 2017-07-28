## Errors about contigs in BAM or VCF files not being properly ordered or sorted

http://gatkforums.broadinstitute.org/gatk/discussion/1328/errors-about-contigs-in-bam-or-vcf-files-not-being-properly-ordered-or-sorted

<p><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/17/14e877060308e4811f8a02c1ca5c85.png" align="right" height="300" /> This is not as common as the &quot;wrong reference build&quot; problem, but it still pops up every now and then: a collaborator gives you a BAM or VCF file that's derived from the correct reference, but for whatever reason the contigs are not sorted in the same order. The GATK can be particular about the <a href="http://www.broadinstitute.org/gatk/guide/article?id=1204">ordering BAM and VCF files</a> so it will fail with an error in this case. </p>
<p>So what do you do? </p>
<hr />
<h3>For BAM files</h3>
<p>You run Picard's <a href="http://broadinstitute.github.io/picard/command-line-overview.html#ReorderSam">ReorderSam</a> tool on your BAM file, using the reference genome dictionary as a template, like this:</p>
<pre><code class="pre_md">java -jar picard.jar ReorderSam \
    I=original.bam \
    O=reordered.bam \
    R=reference.fasta \
    CREATE_INDEX=TRUE</code class="pre_md"></pre>
<p>Where <code>reference.fasta</code> is your genome reference, which <em>must</em> be accompanied by a valid <code>*.dict</code> dictionary file. The <code>CREATE_INDEX</code> argument is optional but useful if you plan to use the resulting file directly with GATK (otherwise you'll need to run another tool to create an index). </p>
<p>Be aware that this tool will drop reads that don't have equivalent contigs in the new reference (potentially bad or not, depending on what you want). If contigs have the same name in the BAM and the new reference, this tool assumes that the alignment of the read in the new BAM is the same. <strong>This is not a liftover tool!</strong></p>
<hr />
<h3>For VCF files</h3>
<p>You run Picard's <a href="http://broadinstitute.github.io/picard/command-line-overview.html#SortVcf">SortVcf</a> tool on your VCF file, using the reference genome dictionary as a template, like this:</p>
<pre><code class="pre_md">java -jar picard.jar SortVcf \
    I=original.vcf \
    O=sorted.vcf \
    SEQUENCE_DICTIONARY=reference.dict </code class="pre_md"></pre>
<p>Where <code>reference.dict</code> is the sequence dictionary of your genome reference.</p>
<p>Note that you may need to delete the index file that gets created automatically for your new VCF by the Picard tool. GATK will automatically regenerate an index file for your VCF. </p>
<h4>Version-specific alert for GATK 3.5</h4>
<p>In version 3.5, we added some beefed-up VCF sequence dictionary validation. Unfortunately, as a side effect of the additional checks, some users have experienced an error that starts with &quot;ERROR MESSAGE: Lexicographically sorted human genome sequence detected in variant.&quot; that is due to unintentional activation of a check that is not necessary. This will be fixed in the next release; in the meantime -U ALLOW_SEQ_DICT_INCOMPATIBILITY can be used (with caution) to override the check.</p>