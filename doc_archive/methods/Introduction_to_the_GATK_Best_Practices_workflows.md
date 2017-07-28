## Introduction to the GATK Best Practices workflows

http://gatkforums.broadinstitute.org/gatk/discussion/4066/introduction-to-the-gatk-best-practices-workflows

<notice><b>This article is part of the Best Practices documentation. See http://www.broadinstitute.org/gatk/guide/best-practices for the full documentation set.</b></notice>
<p>The &quot;GATK Best Practices&quot; are workflow descriptions that provide step-by-step recommendations for getting the best analysis results possible out of high-throughput sequencing data. At present, we provide the following Best Practice workflows:</p>
<ul>
<li><a href="https://www.broadinstitute.org/gatk/guide/best-practices?bpm=DNAseq">Variant Discovery in DNAseq</a></li>
<li><a href="https://www.broadinstitute.org/gatk/guide/best-practices?bpm=RNAseq">Variant Discovery in RNAseq</a></li>
</ul>
<p>These recommendations have been developed by the <a href="http://www.broadinstitute.org/gatk/about/who-we-are">GATK development team</a> over years of analysis work on many of the Broad Institute's sequencing projects, and are applied in the Broad's production pipelines. As a general rule, the command-line arguments and parameters given in the documentation examples are meant to be broadly applicable.</p>
<hr />
<h4>Important notes on context and caveats</h4>
<p>Our testing focuses largely on data from human whole-genome or whole-exome samples sequenced with Illumina technology, so if you are working with different types of data or experimental designs, you may need to adapt certain branches of the workflow, as well as certain parameter selections and values. Unfortunately we are not able to provide official recommendations on how to deal with very different experimental designs or divergent datatypes (such as Ion Torrent).</p>
<p>In addition, the illustrations and tutorials provided in these pages tend to assume a simple experimental design where each sample is used to produce one DNA library that is sequenced separately on one lane of the machine. See the Guide for help dealing with other experimental designs.</p>
<p>Finally, please be aware that several key steps in the Best Practices workflow make use of existing resources such as known variants, which are readily available for humans (we provide several useful resource datasets for download from our FTP server). If no such resources are available for your organism, you may need to bootstrap your own or use alternative methods. We have documented useful methods to do this wherever possible, but be aware than some issues are currently still without a good solution.</p>
<hr />
<notice><b>Important note on GATK versions</b></notice>
<version>
The <a href='http://www.broadinstitute.org/gatk/guide/best-practices'>Best Practices</a> have been updated for GATK version 3. If you are running an older version, you should seriously consider upgrading. For more details about what has changed in each version, please see the <a href='http://www.broadinstitute.org/gatk/guide/version-history'>Version History</a> section. If you cannot upgrade your version of GATK for any reason, please look up the corresponding version of the GuideBook PDF (also in the <a href='http://www.broadinstitute.org/gatk/guide/version-history'>Version History</a> section) to ensure that you are using the appropriate recommendations for your version.</version>