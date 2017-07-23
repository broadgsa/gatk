## Combining variants from different files into one

http://gatkforums.broadinstitute.org/gatk/discussion/53/combining-variants-from-different-files-into-one

<h3>Solutions for combining variant callsets depending on purpose</h3>
<p>There are three main reasons why you might want to combine variants from different files into one, and the tool to use depends on what you are trying to achieve.</p>
<ol>
<li>
<p>The most common case is when you have been parallelizing your variant calling analyses, e.g. running HaplotypeCaller per-chromosome, producing separate VCF files (or gVCF files) per-chromosome. For that case, you can use a tool called <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_CatVariants.php">CatVariants</a> to concatenate the files. There are a few important requirements (e.g. the files should contain all the same samples, and distinct intervals) which you can read about on the tool's documentation page. </p>
</li>
<li>
<p>The second case is when you have been using HaplotypeCaller in <code>-ERC GVCF</code> or <code>-ERC BP_RESOLUTION</code> to call variants on a large cohort, producing many gVCF files. We recommend combining the output gVCF in batches of e.g. 200 before putting them through joint genotyping with GenotypeGVCFs  (for performance reasons), which you can do using <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_CombineGVCFs.php">CombineGVCFs</a>, which is specific for handling gVCF files.</p>
</li>
<li>The third case is when you want to combine variant calls that were produced from the same samples but using different methods, for comparison. For example, if you're evaluating variant calls produced by different variant callers, different workflows, or the same but using different parameters. This produces separate callsets for the same samples, which are then easier to compare if you combine them into a single file. For that purpose, you can use <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_CombineVariants.php">CombineVariants</a>, which is capable of merging VCF records intelligently, treating the same samples as separate or not as desired, combining annotations as appropriate. This is the case that requires the most preparation and forethought because there are many options that may be used to adapt the behavior of the tool.</li>
</ol>
<p>There is also one reason you might want to combine variants from different files into one, that we do not recommend following. That is, if you have produced variant calls from various samples separately, and want to combine them for analysis. This is how people used to do variant analysis on large numbers of samples, but we don't recommend proceeding this way because that workflow suffers from serious methodological flaws. Instead, you should follow our recommendations as laid out in the <a href="https://www.broadinstitute.org/gatk/guide/best-practices">Best Practices</a> documentation.   </p>
<hr />
<h3>Merging records across VCFs with <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_CombineVariants.php">CombineVariants</a></h3>
<p>Here we provide some more information and a worked out example to illustrate the third case because it is less straightforward than the other two.</p>
<p>A key point to understand is that CombineVariants will include a record at every site in all of your input VCF files, and annotate in which input callsets the record is present, pass, or filtered in in the set attribute in the <code>INFO</code> field (see below). In effect, CombineVariants always produces a union of the input VCFs. Any part of the Venn of the N merged VCFs can then be extracted specifically using JEXL expressions on the set attribute using SelectVariants. If you want to extract just the records in common between two VCFs, you would first CombineVariants the two files into a single VCF, and then run SelectVariants to extract the common records with <code>-select 'set == "Intersection"'</code>, as worked out in the detailed example below.</p>
<h4>Handling PASS/FAIL records at the same site in multiple input files</h4>
<p>The <code>-filteredRecordsMergeType</code> argument determines how CombineVariants handles sites where a record is present in multiple VCFs, but it is filtered in some and unfiltered in others, as described in the tool documentation page linked above.</p>
<h4>Understanding the set attribute</h4>
<p>The set property of the <code>INFO</code> field indicates which call set the variant was found in. It can take on a variety of values indicating the exact nature of the overlap between the call sets. Note that the values are generalized for multi-way combinations, but here we describe only the values for 2 call sets being combined.</p>
<ul>
<li>
<p><code>set=Intersection</code> : occurred in both call sets, not filtered out</p>
</li>
<li>
<p><code>set=NAME</code> : occurred in the call set <code>NAME</code> only</p>
</li>
<li>
<p><code>set=NAME1-filteredInNAME</code> : occurred in both call sets, but was not filtered in <code>NAME1</code> but was filtered in <code>NAME2</code></p>
</li>
<li><code>set=filteredInAll</code> : occurred in both call sets, but was filtered out of both</li>
</ul>
<p>For three or more call sets combinations, you can see records like <code>NAME1-NAME2</code> indicating a variant occurred in both <code>NAME1</code> and <code>NAME2</code> but not all sets.</p>
<p>You specify the <code>NAME</code> of a callset is by using the following syntax in your command line: <code>-V:omni 1000G_omni2.5.b37.sites.vcf</code>.</p>
<h4>Emitting minimal VCF output</h4>
<p>You can add the <code>-minimalVCF</code> argument to CombineVariants if you want to eliminate unnecessary information from the <code>INFO</code> field and genotypes. In that case, the only fields emitted will be <code>GT:GQ</code> for genotypes and the <code>keySet</code> for <code>INFO</code>.</p>
<p>An even more extreme output format is <code>-sites_only</code> (a general engine capability listed in the <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_engine_CommandLineGATK.php">CommandLineGATK documentation</a>) where the genotypes for all samples are completely stripped away from the output format.  Enabling this option results in a significant performance speedup as well.</p>
<h4>Requiring sites to be present in a minimum number of callsets</h4>
<p>Sometimes you may want to combine several data sets but you only keep sites that are present in at least 2 of them. To do so, simply add the <code>-minN</code> (or <code>--minimumN</code>) command, followed by an integer if you want to only output records present in at least N input files. In our example, you would add <code>-minN 2</code> to the command line.</p>
<h4>Example: intersecting two VCFs</h4>
<p>In the following example, we use CombineVariants and SelectVariants to obtain only the sites in common between the OMNI 2.5M and HapMap3 sites in the GSA bundle.</p>
<pre><code class="pre_md"># combine the data
java -Xmx2g -jar dist/GenomeAnalysisTK.jar -T CombineVariants -R bundle/b37/human_g1k_v37.fasta -L 1:1-1,000,000 -V:omni bundle/b37/1000G_omni2.5.b37.sites.vcf -V:hm3 bundle/b37/hapmap_3.3.b37.sites.vcf -o union.vcf

# select the intersection
java -Xmx2g -jar dist/GenomeAnalysisTK.jar -T SelectVariants -R ~/Desktop/broadLocal/localData/human_g1k_v37.fasta -L 1:1-1,000,000 -V:variant union.vcf -select 'set == "Intersection";' -o intersect.vcf</code class="pre_md"></pre>
<p>This results in two vcf files, which look like:</p>
<pre><code class="pre_md"># contents of union.vcf
1       990839  SNP1-980702     C       T       .       PASS    AC=150;AF=0.05384;AN=2786;CR=100.0;GentrainScore=0.7267;HW=0.0027632264;set=Intersection
1       990882  SNP1-980745     C       T       .       PASS    CR=99.79873;GentrainScore=0.7403;HW=0.005225421;set=omni
1       990984  SNP1-980847     G       A       .       PASS    CR=99.76005;GentrainScore=0.8406;HW=0.26163524;set=omni
1       992265  SNP1-982128     C       T       .       PASS    CR=100.0;GentrainScore=0.7412;HW=0.0025895447;set=omni
1       992819  SNP1-982682     G       A       .       id50    CR=99.72961;GentrainScore=0.8505;HW=4.811053E-17;set=FilteredInAll
1       993987  SNP1-983850     T       C       .       PASS    CR=99.85935;GentrainScore=0.8336;HW=9.959717E-28;set=omni
1       994391  rs2488991       G       T       .       PASS    AC=1936;AF=0.69341;AN=2792;CR=99.89378;GentrainScore=0.7330;HW=1.1741E-41;set=filterInomni-hm3
1       996184  SNP1-986047     G       A       .       PASS    CR=99.932205;GentrainScore=0.8216;HW=3.8830226E-6;set=omni
1       998395  rs7526076       A       G       .       PASS    AC=2234;AF=0.80187;AN=2786;CR=100.0;GentrainScore=0.8758;HW=0.67373306;set=Intersection
1       999649  SNP1-989512     G       A       .       PASS    CR=99.93262;GentrainScore=0.7965;HW=4.9767335E-4;set=omni

# contents of intersect.vcf
1       950243  SNP1-940106     A       C       .       PASS    AC=826;AF=0.29993;AN=2754;CR=97.341675;GentrainScore=0.7311;HW=0.15148845;set=Intersection
1       957640  rs6657048       C       T       .       PASS    AC=127;AF=0.04552;AN=2790;CR=99.86667;GentrainScore=0.6806;HW=2.286109E-4;set=Intersection
1       959842  rs2710888       C       T       .       PASS    AC=654;AF=0.23559;AN=2776;CR=99.849;GentrainScore=0.8072;HW=0.17526293;set=Intersection
1       977780  rs2710875       C       T       .       PASS    AC=1989;AF=0.71341;AN=2788;CR=99.89077;GentrainScore=0.7875;HW=2.9912625E-32;set=Intersection
1       985900  SNP1-975763     C       T       .       PASS    AC=182;AF=0.06528;AN=2788;CR=99.79926;GentrainScore=0.8374;HW=0.017794203;set=Intersection
1       987200  SNP1-977063     C       T       .       PASS    AC=1956;AF=0.70007;AN=2794;CR=99.45917;GentrainScore=0.7914;HW=1.413E-42;set=Intersection
1       987670  SNP1-977533     T       G       .       PASS    AC=2485;AF=0.89196;AN=2786;CR=99.51427;GentrainScore=0.7005;HW=0.24214932;set=Intersection
1       990417  rs2465136       T       C       .       PASS    AC=1113;AF=0.40007;AN=2782;CR=99.7599;GentrainScore=0.8750;HW=8.595538E-5;set=Intersection
1       990839  SNP1-980702     C       T       .       PASS    AC=150;AF=0.05384;AN=2786;CR=100.0;GentrainScore=0.7267;HW=0.0027632264;set=Intersection
1       998395  rs7526076       A       G       .       PASS    AC=2234;AF=0.80187;AN=2786;CR=100.0;GentrainScore=0.8758;HW=0.67373306;set=Intersection</code class="pre_md"></pre>