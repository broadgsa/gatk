## How should I cite GATK in my own publications?

http://gatkforums.broadinstitute.org/gatk/discussion/6201/how-should-i-cite-gatk-in-my-own-publications

<p>To date we have published three papers on GATK (citation details below). The ideal way to cite the GATK is to use all as a triple citation, as in:</p>
<blockquote>
<p>We sequenced 10 samples on 10 lanes on an Illumina HiSeq 2000, aligned the resulting reads to the hg19 reference genome with BWA (Li &amp; Durbin), applied GATK <strong>(McKenna <em>et al.</em>, 2010)</strong> base quality score recalibration, indel realignment, duplicate removal, and performed SNP and INDEL discovery and genotyping across all 10 samples simultaneously using standard hard filtering parameters or variant quality score recalibration according to GATK Best Practices recommendations <strong>(DePristo <em>et al.</em>, 2011; Van der Auwera <em>et al.</em>, 2013)</strong>.</p>
</blockquote>
<hr />
<h3>McKenna <em>et al.</em> 2010 : Original description of the GATK framework</h3>
<p>The first GATK paper covers the computational philosophy underlying the GATK and is a good citation for the GATK in general.</p>
<p><strong>The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data</strong> McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA, 2010 <em>GENOME RESEARCH 20:1297-303</em> </p>
<p><a href="http://dx.doi.org/10.1101/gr.107524.110">Article</a> | <a href="http://www.ncbi.nlm.nih.gov/pubmed?term=20644199">Pubmed</a></p>
<hr />
<h3>DePristo <em>et al.</em> 2011 : First incarnation of the Best Practices workflow</h3>
<p>The second GATK paper describes in more detail some of the key tools commonly used in the GATK for high-throughput sequencing data processing and variant discovery. The paper covers base quality score recalibration, indel realignment, SNP calling with UnifiedGenotyper, variant quality score recalibration and their application to deep whole genome, whole exome, and low-pass multi-sample calling. This is a good citation if you use the GATK for variant discovery. </p>
<p><strong>A framework for variation discovery and genotyping using next-generation DNA sequencing data</strong> DePristo M, Banks E, Poplin R, Garimella K, Maguire J, Hartl C, Philippakis A, del Angel G, Rivas MA, Hanna M, McKenna A, Fennell T, Kernytsky A, Sivachenko A, Cibulskis K, Gabriel S, Altshuler D, Daly M, 2011 <em>NATURE GENETICS 43:491-498</em> </p>
<p><a href="http://dx.doi.org/10.1038/ng.806">Article</a> | <a href="http://www.ncbi.nlm.nih.gov/pubmed?term=21478889">Pubmed</a></p>
<p>Note that the workflow described in this paper corresponds to the version 1.x to 2.x best practices. Some key steps for variant discovery have been significantly modified in later versions (3.x onwards). This paper should not be used as a definitive guide to variant discovery with GATK. For that, please see our online documentation guide.</p>
<hr />
<h3>Van der Auwera <em>et al.</em> 2013 : Hands-on tutorial with step-by-step explanations</h3>
<p>The third GATK paper describes the Best Practices for Variant Discovery (version 2.x). It is intended mainly as a learning resource for first-time users and as a protocol reference. This is a good citation to include in a Materials and Methods section. </p>
<p><strong>From FastQ Data to High-Confidence Variant Calls: The Genome Analysis Toolkit Best Practices Pipeline</strong> Van der Auwera GA, Carneiro M, Hartl C, Poplin R, del Angel G, Levy-Moonshine A, Jordan T, Shakir K, Roazen D, Thibault J, Banks E, Garimella K, Altshuler D, Gabriel S, DePristo M, 2013 <em>CURRENT PROTOCOLS IN BIOINFORMATICS 43:11.10.1-11.10.33</em> </p>
<p><a href="http://dx.doi.org/10.1002/0471250953.bi1110s43">Article</a> | <a href="http://www.ncbi.nlm.nih.gov/pubmed/?term=25431634">PubMed</a></p>
<p>Remember that as our work continues and our Best Practices recommendations evolve, specific command lines, argument values and even tool choices described in the paper become obsolete. Be sure to always refer to our Best Practices documentation for the most up-to-date and version-appropriate recommendations.</p>