## Which datasets should I use for reviewing or benchmarking purposes?

http://gatkforums.broadinstitute.org/gatk/discussion/1292/which-datasets-should-i-use-for-reviewing-or-benchmarking-purposes

<h2>New WGS and WEx CEU trio BAM files</h2>
<p>We have sequenced at the Broad Institute and released to the 1000 Genomes Project the following datasets for the three members of the CEU trio (NA12878, NA12891 and NA12892):</p>
<ul>
<li>WEx (150x) sequence</li>
<li>WGS (>60x) sequence </li>
</ul>
<p>This is better data to work with than the original DePristo et al. BAMs files, so we recommend you download and analyze these files if you are looking for complete, large-scale data sets to evaluate the GATK or other tools.  </p>
<p>Here's the rough library properties of the BAMs:</p>
<p><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/bc/03c8ec93b4e4f16f14cd9f4d1ae77f.jpeg" alt="CEU trio BAM libraries" /></p>
<p>These data files can be downloaded from the <a href="ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20120117_ceu_trio_b37_decoy/">1000 Genomes DCC</a></p>
<h2>NA12878 Datasets from DePristo et al. (2011) Nature Genetics</h2>
<p>Here are the datasets we used in the GATK paper cited below.</p>
<p><strong>DePristo M, Banks E, Poplin R, Garimella K, Maguire J, Hartl C, Philippakis A, del Angel G, Rivas MA, Hanna M, McKenna A, Fennell T, Kernytsky A, Sivachenko A, Cibulskis K, Gabriel S, Altshuler D and Daly, M (2011).</strong> A framework for variation discovery and genotyping using next-generation DNA sequencing data. <em>Nature Genetics.</em> <strong>43:</strong>491-498.</p>
<p>Some of the BAM and VCF files are currently hosted by the NCBI:
ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20101201_cg_NA12878/</p>
<ul>
<li>NA12878.hiseq.wgs.bwa.recal.bam -- BAM file for NA12878 HiSeq whole genome</li>
<li>NA12878.hiseq.wgs.bwa.raw.bam Raw reads (in BAM format, see below)</li>
<li>NA12878.ga2.exome.maq.recal.bam -- BAM file for NA12878 GenomeAnalyzer II whole exome (hg18)</li>
<li>NA12878.ga2.exome.maq.raw.bam Raw reads (in BAM format, see below)</li>
<li>NA12878.hiseq.wgs.vcf.gz -- SNP calls for NA12878 HiSeq whole genome (hg18)</li>
<li>NA12878.ga2.exome.vcf.gz -- SNP calls for NA12878 GenomeAnalyzer II whole exome (hg18)</li>
<li>BAM files for CEU + NA12878 whole genome (b36). These are the standard BAM files for the 1000 Genomes pilot CEU samples plus a 4x downsampled version of NA12878 from the pilot 2 data set, available in the DePristoNatGenet2011 directory of the <a href="http://www.broadinstitute.org/gatk/guide/article?id=1215">GSA FTP Server</a></li>
<li>SNP calls for CEU + NA12878 whole genome (b36) are available in the DePristoNatGenet2011 directory of the <a href="http://www.broadinstitute.org/gatk/guide/article?id=1215">GSA FTP Server</a></li>
<li>Crossbow comparison SNP calls are available in the DePristoNatGenet2011 directory of the <a href="http://www.broadinstitute.org/gatk/guide/article?id=1215">GSA FTP Server</a> as crossbow.filtered.vcf. The raw calls can be viewed by ignoring the <code>FILTER</code> field status</li>
<li>whole_exome_agilent_designed_120.Homo_sapiens_assembly18.targets.interval_list <code>-- targets</code> used in the analysis of the exome capture data</li>
</ul>
<p>Please note that we have not collected the indel calls for the paper, as these are only used for filtering SNPs near indels. If you want to call accurate indels, please use the new GATK indel caller in the Unified Genotyper.</p>
<h3>Warnings</h3>
<p>Both the GATK and the sequencing technologies have improved significantly since the analyses performed in this paper.</p>
<ul>
<li>
<p>If you are conducting a review today, we would recommend that the newest version of the GATK, which performs much better than the version described in the paper. Moreover, we would also recommend one use the newest version of Crossbow as well, in case they have improved things. The GATK calls for NA12878 from the paper (above) will give one a good idea what a good call set looks like whole-genome or whole-exome.</p>
</li>
<li>The data sets used in the paper are no longer state-of-the-art. The WEx BAM is GAII data aligned with MAQ on hg18, but a state-of-the-art data set would use HiSeq and BWA on hg19. Even the 64x HiSeq WG data set is already more than one year old. For a better assessment, we would recommend you use a newer data set for these samples, if you have the capacity to generate it. This applies less to the WG NA12878 data, which is pretty good, but the NA12878 WEx from the paper is nearly 2 years old now and notably worse than our most recent data sets.</li>
</ul>
<p>Obviously, this was an annoyance for us as well, as it would have been nice to use a state-of-the-art data set for the WEx. But we decided to freeze the data used for analysis to actually finish this paper.</p>
<h3>How do I get the raw FASTQ file from a BAM?</h3>
<p>If you want the raw, machine output for the data analyzed in the GATK framework paper, obtain the raw BAM files above and convert them from SAM to FASTQ using the Picard tool <a href="http://picard.sourceforge.net/command-line-overview.shtml#SamToFastq">SamToFastq</a>.</p>