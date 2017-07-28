## Adding Genomic Annotations Using SnpEff and VariantAnnotator

http://gatkforums.broadinstitute.org/gatk/discussion/50/adding-genomic-annotations-using-snpeff-and-variantannotator

<h3>This article is out of date and no longer applicable. At this time, we do not provide support for performing functional annotation. Programs that we are aware of and that our collaborators use successfully include Oncotator and Variant Effect Predictor (VEP).</h3>
<hr />
<p><em>Our testing has shown that not all combinations of snpEff/database versions produce high-quality results. Be sure to read this document completely to familiarize yourself with our recommended best practices BEFORE running snpEff.</em></p>
<h3>Introduction</h3>
<p>Until recently we were using an in-house annotation tool for genomic annotation, but the burden of keeping the database current and our lack of ability to annotate indels has led us to employ the use of a third-party tool instead. After reviewing many external tools (including annoVar, VAT, and Oncotator), we decided that <a href="http://snpeff.sourceforge.net/">SnpEff</a> best meets our needs as it accepts VCF files as input, can annotate a full exome callset (including indels) in seconds, and provides continually-updated transcript databases. We have implemented support in the GATK for parsing the output from the SnpEff tool and annotating VCFs with the information provided in it. </p>
<h3>SnpEff Setup and Usage</h3>
<p>Download the SnpEff core program. If you want to be able to run VariantAnnotator on the SnpEff output, you'll need to download a version of SnpEff that VariantAnnotator supports from <a href="http://sourceforge.net/projects/snpeff/files/">this page</a> (currently supported versions are listed below). If you just want the most recent version of SnpEff and don't plan to run VariantAnnotator on its output, you can get it from <a href="http://snpeff.sourceforge.net/download.html">here</a>.</p>
<p>After unzipping the core program, open the file snpEff.config in a text editor, and change the &quot;database_repository&quot; line to the following:</p>
<pre><code class="pre_md">database_repository = http://sourceforge.net/projects/snpeff/files/databases/</code class="pre_md"></pre>
<p>Then, download one or more databases using SnpEff's built-in download command:</p>
<pre><code class="pre_md">java -jar snpEff.jar download GRCh37.64</code class="pre_md"></pre>
<p>You can find a list of available databases <a href="http://snpeff.sourceforge.net/download.html">here</a>. The human genome databases have <strong>GRCh</strong> or <strong>hg</strong> in their names. You can also download the databases directly from the SnpEff website, if you prefer.</p>
<p>The download command by default puts the databases into a subdirectory called <strong>data</strong> within the directory containing the SnpEff jar file. If you want the databases in a different directory, you'll need to edit the <code>data_dir</code> entry in the file <code>snpEff.config</code> to point to the correct directory.</p>
<p>Run SnpEff on the file containing your variants, and redirect its output to a file. SnpEff supports many input file formats including VCF 4.1, BED, and SAM pileup. Full details and command-line options can be found on the <a href="http://snpeff.sourceforge.net/">SnpEff home page</a>.</p>
<h3>Supported SnpEff Versions</h3>
<p>If you want to take advantage of SnpEff integration in the GATK, you'll need to run SnpEff version *<em>2.0.5</em>. <em>Note: newer versions are currently unsupported by the GATK, as we haven't yet had the reources to test it.</em></p>
<h3>Current Recommended Best Practices When Running SnpEff</h3>
<p>These best practices are based on our analysis of various snpEff/database versions as described in detail in the <strong>Analysis of SnpEff Annotations Across Versions</strong> section below.</p>
<ul>
<li>
<p>We recommend using only the <strong>GRCh37.64</strong> database with SnpEff 2.0.5. The more recent GRCh37.65 database produces many false-positive Missense annotations due to a regression in the ENSEMBL Release 65 GTF file used to build the database. This regression has been acknowledged by ENSEMBL and is supposedly fixed as of 1-30-2012; however as we have not yet tested the fixed version of the database we continue to recommend using only GRCh37.64 for now.</p>
</li>
<li>
<p>We recommend always running with <code>-onlyCoding true</code> with human databases (eg., the GRCh37.<em> databases). Setting <code>-onlyCoding false</code> causes snpEff to report all transcripts as if they were coding (even if they're not), which can lead to nonsensical results. The <code>-onlyCoding false</code> option should </em>only* be used with databases that lack protein coding information.</p>
</li>
<li>Do not trust annotations from versions of snpEff prior to 2.0.4. Older versions of snpEff (such as 2.0.2) produced many incorrect annotations due to the presence of a certain number of nonsensical transcripts in the underlying ENSEMBL databases. Newer versions of snpEff filter out such transcripts.</li>
</ul>
<h3>Analyses of SnpEff Annotations Across Versions</h3>
<p>See our analysis of the SNP annotations produced by snpEff across various snpEff/database versions <a href="http://www.broadinstitute.org/gatk/media/docs/SnpEff_snps_comparison_of_available_versions.pdf">here</a>.</p>
<ul>
<li>
<p>Both snpEff 2.0.2 + GRCh37.63 and snpEff 2.0.5 + GRCh37.65 produce an abnormally high Missense:Silent ratio, with elevated levels of Missense mutations across the entire spectrum of allele counts. They also have a relatively low (~70%) level of concordance with the 1000G Gencode annotations when it comes to Silent mutations. This suggests that these combinations of snpEff/database versions incorrectly annotate many Silent mutations as Missense.</p>
</li>
<li>snpEff 2.0.4 RC3 + GRCh37.64 and snpEff 2.0.5 + GRCh37.64 produce a Missense:Silent ratio in line with expectations, and have a very high (~97%-99%) level of concordance with the 1000G Gencode annotations across all categories.</li>
</ul>
<p>See our comparison of SNP annotations produced using the GRCh37.64 and GRCh37.65 databases with snpEff 2.0.5 <a href="http://www.broadinstitute.org/gatk/media/docs/SnpEff_snps_ensembl_64_vs_65.pdf">here</a></p>
<ul>
<li>
<p>The GRCh37.64 database gives good results on the condition that you run snpEff with the <code>-onlyCoding true</code> option. The <code>-onlyCoding false</code> option causes snpEff to mark <em>all</em> transcripts as coding, and so produces many false-positive Missense annotations.</p>
</li>
<li>The GRCh37.65 database gives results that are as poor as those you get with the <code>-onlyCoding false</code> option on the GRCh37.64 database. This is due to a regression in the ENSEMBL release 65 GTF file used to build snpEff's GRCh37.65 database. The regression has been acknowledged by ENSEMBL and is due to be fixed shortly.</li>
</ul>
<p>See our analysis of the INDEL annotations produced by snpEff across snpEff/database versions <a href="http://www.broadinstitute.org/gatk/media/docs/SnpEff_indels.pdf">here</a></p>
<ul>
<li>snpEff's indel annotations are highly concordant with those of a high-quality set of genomic annotations from the 1000 Genomes project. This is true across all snpEff/database versions tested.</li>
</ul>
<h3>Example SnpEff Usage with a VCF Input File</h3>
<p>Below is an example of how to run SnpEff version 2.0.5 with a VCF input file and have it write its output in VCF format as well. Notice that you need to explicitly specify the database you want to use (in this case, GRCh37.64). This database must be present in a directory of the same name within the <code>data_dir</code> as defined in <code>snpEff.config</code>.</p>
<pre><code class="pre_md">java -Xmx4G -jar snpEff.jar eff -v -onlyCoding true -i vcf -o vcf GRCh37.64 1000G.exomes.vcf &gt; snpEff_output.vcf</code class="pre_md"></pre>
<p>In this mode, SnpEff aggregates all effects associated with each variant record together into a single INFO field annotation with the key EFF. The general format is:</p>
<pre><code class="pre_md">EFF=Effect1(Information about Effect1),Effect2(Information about Effect2),etc.</code class="pre_md"></pre>
<p>And here is the precise layout with all the subfields:</p>
<pre><code class="pre_md">EFF=Effect1(Effect_Impact|Effect_Functional_Class|Codon_Change|Amino_Acid_Change|Gene_Name|Gene_BioType|Coding|Transcript_ID|Exon_ID),Effect2(etc...</code class="pre_md"></pre>
<p>It's also possible to get SnpEff to output in a (non-VCF) text format with one Effect per line. See the <a href="http://snpeff.sourceforge.net/">SnpEff home page</a> for full details.</p>
<h3>Adding SnpEff Annotations using VariantAnnotator</h3>
<p>Once you have a SnpEff output VCF file, you can use the VariantAnnotator walker to add SnpEff annotations based on that output to the input file you ran SnpEff on.</p>
<p>There are two different options for doing this:</p>
<h4>Option 1: Annotate with only the highest-impact effect for each variant</h4>
<p><em>NOTE: This option works only with supported SnpEff versions as explained above. VariantAnnotator run as described below will refuse to parse SnpEff output files produced by other versions of the tool, or which lack a SnpEff version number in their header.</em></p>
<p>The default behavior when you run VariantAnnotator on a SnpEff output file is to parse the complete set of effects resulting from the current variant, select the most biologically-significant effect, and add annotations for just that effect to the INFO field of the VCF record for the current variant. This is the mode we plan to use in our Production Data-Processing Pipeline.</p>
<p>When selecting the most biologically-significant effect associated with the current variant, VariantAnnotator does the following:</p>
<ul>
<li>
<p>Prioritizes the effects according to the categories (in order of decreasing precedence) &quot;High-Impact&quot;, &quot;Moderate-Impact&quot;, &quot;Low-Impact&quot;, and &quot;Modifier&quot;, and always selects one of the effects from the highest-priority category. For example, if there are three moderate-impact effects and two high-impact effects resulting from the current variant, the annotator will choose one of the high-impact effects and add annotations based on it. See below for a full list of the effects arranged by category.</p>
</li>
<li>
<p>Within each category, ties are broken using the functional class of each effect (in order of precedence: NONSENSE, MISSENSE, SILENT, or NONE). For example, if there is both a NON_SYNONYMOUS_CODING (MODERATE-impact, MISSENSE) and a CODON_CHANGE (MODERATE-impact, NONE) effect associated with the current variant, the annotator will select the NON_SYNONYMOUS_CODING effect. This is to allow for more accurate counts of the total number of sites with NONSENSE/MISSENSE/SILENT mutations. See below for a description of the functional classes SnpEff associates with the various effects.</p>
</li>
<li>Effects that are within a non-coding region are always considered lower-impact than effects that are within a coding region.</li>
</ul>
<p>Example Usage:</p>
<pre><code class="pre_md">java -jar dist/GenomeAnalysisTK.jar \
     -T VariantAnnotator \
     -R /humgen/1kg/reference/human_g1k_v37.fasta \
     -A SnpEff \       
     --variant 1000G.exomes.vcf \        (file to annotate)
     --snpEffFile snpEff_output.vcf \    (SnpEff VCF output file generated by running SnpEff on the file to annotate)
     -L 1000G.exomes.vcf \
     -o out.vcf</code class="pre_md"></pre>
<p>VariantAnnotator adds some or all of the following INFO field annotations to each variant record:</p>
<ul>
<li><code>SNPEFF_EFFECT</code> - The highest-impact effect resulting from the current variant (or one of the highest-impact effects, if there is a tie)</li>
<li><code>SNPEFF_IMPACT</code> - Impact of the highest-impact effect resulting from the current variant (<code>HIGH</code>, <code>MODERATE</code>, <code>LOW</code>, or <code>MODIFIER</code>)</li>
<li><code>SNPEFF_FUNCTIONAL_CLASS</code> - Functional class of the highest-impact effect resulting from the current variant (<code>NONE</code>, <code>SILENT</code>, <code>MISSENSE</code>, or <code>NONSENSE</code>)</li>
<li><code>SNPEFF_CODON_CHANGE</code> - Old/New codon for the highest-impact effect resulting from the current variant</li>
<li><code>SNPEFF_AMINO_ACID_CHANGE</code> - Old/New amino acid for the highest-impact effect resulting from the current variant</li>
<li><code>SNPEFF_GENE_NAME</code> - Gene name for the highest-impact effect resulting from the current variant</li>
<li><code>SNPEFF_GENE_BIOTYPE</code> - Gene biotype for the highest-impact effect resulting from the current variant</li>
<li><code>SNPEFF_TRANSCRIPT_ID</code> - Transcript ID for the highest-impact effect resulting from the current variant</li>
<li><code>SNPEFF_EXON_ID</code> - Exon ID for the highest-impact effect resulting from the current variant</li>
</ul>
<p>Example VCF records annotated using SnpEff and VariantAnnotator:</p>
<pre><code class="pre_md">1   874779  .   C   T   279.94  . AC=1;AF=0.0032;AN=310;BaseQRankSum=-1.800;DP=3371;Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=1.4493;InbreedingCoeff=-0.0045;
MQ=54.49;MQ0=10;MQRankSum=0.982;QD=13.33;ReadPosRankSum=-0.060;SB=-120.09;SNPEFF_AMINO_ACID_CHANGE=G215;SNPEFF_CODON_CHANGE=ggC/ggT;
SNPEFF_EFFECT=SYNONYMOUS_CODING;SNPEFF_EXON_ID=exon_1_874655_874840;SNPEFF_FUNCTIONAL_CLASS=SILENT;SNPEFF_GENE_BIOTYPE=protein_coding;SNPEFF_GENE_NAME=SAMD11;
SNPEFF_IMPACT=LOW;SNPEFF_TRANSCRIPT_ID=ENST00000342066

1   874816  .   C   CT  2527.52 .   AC=15;AF=0.0484;AN=310;BaseQRankSum=-11.876;DP=4718;FS=48.575;HRun=1;HaplotypeScore=91.9147;InbreedingCoeff=-0.0520;
MQ=53.37;MQ0=6;MQRankSum=-1.388;QD=5.92;ReadPosRankSum=-1.932;SB=-741.06;SNPEFF_EFFECT=FRAME_SHIFT;SNPEFF_EXON_ID=exon_1_874655_874840;
SNPEFF_FUNCTIONAL_CLASS=NONE;SNPEFF_GENE_BIOTYPE=protein_coding;SNPEFF_GENE_NAME=SAMD11;SNPEFF_IMPACT=HIGH;SNPEFF_TRANSCRIPT_ID=ENST00000342066</code class="pre_md"></pre>
<h4>Option 2: Annotate with all effects for each variant</h4>
<p>VariantAnnotator also has the ability to take the EFF field from the SnpEff VCF output file containing all the effects aggregated together and copy it verbatim into the VCF to annotate.</p>
<p>Here's an example of how to do this:</p>
<pre><code class="pre_md">java -jar dist/GenomeAnalysisTK.jar \
     -T VariantAnnotator \
     -R /humgen/1kg/reference/human_g1k_v37.fasta \      
     -E resource.EFF \
     --variant 1000G.exomes.vcf \      (file to annotate)
     --resource snpEff_output.vcf \    (SnpEff VCF output file generated by running SnpEff on the file to annotate)
     -L 1000G.exomes.vcf \
     -o out.vcf</code class="pre_md"></pre>
<p>Of course, in this case you can also use the VCF output by SnpEff directly, but if you are using VariantAnnotator for other purposes anyway the above might be useful.</p>
<h3>List of Genomic Effects</h3>
<p>Below are the possible genomic effects recognized by SnpEff, grouped by biological impact. Full descriptions of each effect are available on <a href="http://snpeff.sourceforge.net/faq.html">this page</a>.</p>
<h4>High-Impact Effects</h4>
<ul>
<li>SPLICE_SITE_ACCEPTOR</li>
<li>SPLICE_SITE_DONOR</li>
<li>START_LOST</li>
<li>EXON_DELETED</li>
<li>FRAME_SHIFT</li>
<li>STOP_GAINED</li>
<li>STOP_LOST</li>
</ul>
<h4>Moderate-Impact Effects</h4>
<ul>
<li>NON_SYNONYMOUS_CODING</li>
<li>CODON_CHANGE <i>(note: this effect is used by SnpEff only for MNPs, not SNPs)</i></li>
<li>CODON_INSERTION</li>
<li>CODON_CHANGE_PLUS_CODON_INSERTION</li>
<li>CODON_DELETION</li>
<li>CODON_CHANGE_PLUS_CODON_DELETION</li>
<li>UTR_5_DELETED</li>
<li>UTR_3_DELETED</li>
</ul>
<h4>Low-Impact Effects</h4>
<ul>
<li>SYNONYMOUS_START</li>
<li>NON_SYNONYMOUS_START</li>
<li>START_GAINED</li>
<li>SYNONYMOUS_CODING</li>
<li>SYNONYMOUS_STOP</li>
<li>NON_SYNONYMOUS_STOP</li>
</ul>
<h4>Modifiers</h4>
<ul>
<li>NONE</li>
<li>CHROMOSOME</li>
<li>CUSTOM</li>
<li>CDS</li>
<li>GENE</li>
<li>TRANSCRIPT</li>
<li>EXON</li>
<li>INTRON_CONSERVED</li>
<li>UTR_5_PRIME</li>
<li>UTR_3_PRIME</li>
<li>DOWNSTREAM</li>
<li>INTRAGENIC</li>
<li>INTERGENIC</li>
<li>INTERGENIC_CONSERVED</li>
<li>UPSTREAM</li>
<li>REGULATION</li>
<li>INTRON</li>
</ul>
<h3>Functional Classes</h3>
<p>SnpEff assigns a functional class to certain effects, in addition to an impact:</p>
<ul>
<li><code>NONSENSE</code>: assigned to point mutations that result in the creation of a new stop codon</li>
<li><code>MISSENSE</code>: assigned to point mutations that result in an amino acid change, but not a new stop codon</li>
<li><code>SILENT</code>: assigned to point mutations that result in a codon change, but not an amino acid change or new stop codon</li>
<li><code>NONE</code>: assigned to all effects that don't fall into any of the above categories (including all events larger than a point mutation)</li>
</ul>
<p>The GATK prioritizes effects with functional classes over effects of equal impact that lack a functional class when selecting the most significant effect in VariantAnnotator. This is to enable accurate counts of NONSENSE/MISSENSE/SILENT sites.</p>