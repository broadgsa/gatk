## Errors about read group (RG) information

http://gatkforums.broadinstitute.org/gatk/discussion/59/errors-about-read-group-rg-information

<h3>What are read groups?</h3>
<p>See the <a href="http://www.broadinstitute.org/gatk/guide/article?id=6472">Dictionary entry</a> on read groups.</p>
<h3>Errors about missing or undefined read groups</h3>
<p>As detailed in the FAQs about input requirements, GATK expects all read groups appearing in the read data to be specified in the file header, and will fail with an error if it does not find that information (whether there is no read group information in the file, or a subset of reads do not have read groups). </p>
<p>Typically you should read group information when you perform the original alignments (with e.g. BWA, which has an option to do so). So what do you do if you forgot to do that, and you don't want to have to rerun BWA all over again? </p>
<h3>Solution</h3>
<p>You can use a Picard tool called <a href="https://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups">AddOrReplaceReadGroups</a> to add the missing information to your input file.</p>
<p>Here's an example:</p>
<pre><code class="pre_md"># throws an error
java -jar GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R reference.fasta \
    -I reads_without_RG.bam \
    -o output.vcf

# fix the read groups
java -jar picard.jar AddOrReplaceReadGroups \
    I= reads_without_RG.bam \
    O=  reads_with_RG.bam \
    SORT_ORDER=coordinate \
    RGID=foo \
    RGLB=bar \
    RGPL=illumina \
    RGSM=Sample1 \
    CREATE_INDEX=True

# runs without error
java -jar GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R reference.fasta \
    -I reads_with_RG.bam \
    -o output.vcf</code class="pre_md"></pre>
<p>Note that if you don't know what information to put in the read groups, you should ask whoever performed the sequencing or provided the BAM to give you the metadata you need.</p>