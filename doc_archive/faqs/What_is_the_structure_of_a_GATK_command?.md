## What is the structure of a GATK command?

http://gatkforums.broadinstitute.org/gatk/discussion/4669/what-is-the-structure-of-a-gatk-command

<h4>Overview</h4>
<p>This document describes how GATK commands are structured and how to add arguments to basic command examples.</p>
<hr />
<h3>Basic java syntax</h3>
<p>Commands for GATK always follow the same basic syntax: </p>
<pre><code class="pre_md">java [Java arguments] -jar GenomeAnalysisTK.jar [GATK arguments]</code class="pre_md"></pre>
<p>The core of the command is <code>java -jar GenomeAnalysisTK.jar</code>, which starts up the GATK program in a Java Virtual Machine (JVM). Any additional java-specific arguments (such as -Xmx to increase memory allocation) should be inserted between <code>java</code> and <code>-jar</code>, like this:</p>
<pre><code class="pre_md">java -Xmx4G -jar GenomeAnalysisTK.jar [GATK arguments]</code class="pre_md"></pre>
<p>The order of arguments between <code>java</code> and <code>-jar</code> is not important.</p>
<hr />
<h3>GATK arguments</h3>
<p>There are two universal arguments that are required for every GATK command (with very few exceptions, the <code>clp</code>-type utilities), <code>-R</code> for Reference (e.g. <code>-R human_b37.fasta</code>) and <code>-T</code> for Tool name (e.g. <code>-T HaplotypeCaller</code>).</p>
<p>Additional arguments fall in two categories: </p>
<ul>
<li>
<p>Engine arguments like <code>-L</code> (for specifying a list of intervals) which can be given to all tools and are technically optional but may be effectively required at certain steps for specific analytical designs (e.g. the <code>-L</code> argument for calling variants on exomes);</p>
</li>
<li>Tool-specific arguments which may be required, like <code>-I</code> (to provide an input file containing sequence reads to tools that process BAM files) or optional, like <code>-alleles</code> (to provide a list of known alleles for genotyping). </li>
</ul>
<p>The ordering of GATK arguments is not important, but we recommend always passing the tool name (<code>-T</code>) and reference (<code>-R</code>) first for consistency. It is also a good idea to consistently order arguments by some kind of logic in order to make it easy to compare different commands over the course of a project. It’s up to you to choose what that logic should be.</p>
<p>All available engine and tool-specific arguments are listed in the <a href="https://www.broadinstitute.org/gatk/guide/tooldocs">tool documentation section</a>. Arguments typically have both a long name (prefixed by <code>--</code>) and a short name (prefixed by <code>-</code>). The GATK command line parser recognizes both equally, so you can use whichever you prefer, depending on whether you prefer commands to be more verbose or more succinct. </p>
<p>Finally, a note about flags. Flags are arguments that have boolean values, i.e. TRUE or FALSE. They are typically used to enable or disable specific features; for example, <code>--keep_program_records</code> will make certain GATK tools output additional information in the BAM header that would be omitted otherwise. In GATK, all flags are set to FALSE by default, so if you want to set one to TRUE, all you need to do is add the flag name to the command. You don't need to specify an actual value.</p>
<hr />
<h3>Examples of complete GATK command lines</h3>
<p>This is a very simple command that runs HaplotypeCaller in default mode on a single input BAM file containing sequence data and outputs a VCF file containing raw variants.</p>
<pre><code class="pre_md">java -Xmx4G -jar GenomeAnalysisTK.jar -R human_b37.fasta -T HaplotypeCaller -I sample1.bam -o raw_variants.vcf</code class="pre_md"></pre>
<p>If the data is from exome sequencing, we should additionally provide the exome targets using the <code>-L</code> argument:</p>
<pre><code class="pre_md">java -Xmx4G -jar GenomeAnalysisTK.jar -R human_b37.fasta -T HaplotypeCaller -I sample1.bam -o raw_variants.vcf -L exome_intervals.list</code class="pre_md"></pre>
<p>If we just want to genotype specific sites of interest using known alleles based on results from a previous study, we can change the HaplotypeCaller’s genotyping mode using <code>-gt_mode</code>, provide those alleles using <code>-alleles</code>, and restrict the analysis to just those sites using <code>-L</code>:</p>
<pre><code class="pre_md">java -Xmx4G -jar GenomeAnalysisTK.jar -R human_b37.fasta -T HaplotypeCaller -I sample1.bam -o raw_variants.vcf -L known_alleles.vcf -alleles known_alleles.vcf -gt_mode GENOTYPE_GIVEN_ALLELES</code class="pre_md"></pre>
<p>For more examples of commands and for specific tool commands, see the <a href="https://www.broadinstitute.org/gatk/guide/tooldocs">tool documentation section</a>.</p>