## What do the VariantEval modules do?

http://gatkforums.broadinstitute.org/gatk/discussion/2361/what-do-the-varianteval-modules-do

<p>VariantEval accepts two types of modules: stratification and evaluation modules.</p>
<ul>
<li>Stratification modules will stratify (group) the variants based on certain properties. </li>
<li>Evaluation modules will compute certain metrics for the variants</li>
</ul>
<h3>CpG</h3>
<p>CpG is a three-state stratification:</p>
<ul>
<li>The locus is a CpG site (&quot;CpG&quot;)</li>
<li>The locus is not a CpG site (&quot;non_CpG&quot;)</li>
<li>The locus is either a CpG or not a CpG site (&quot;all&quot;)</li>
</ul>
<p>A CpG site is defined as a site where the reference base at a locus is a C and the adjacent reference base in the 3' direction is a G.</p>
<h3>EvalRod</h3>
<p>EvalRod is an N-state stratification, where N is the number of eval rods bound to VariantEval.</p>
<h3>Sample</h3>
<p>Sample is an N-state stratification, where N is the number of samples in the eval files.</p>
<h3>Filter</h3>
<p>Filter is a three-state stratification:</p>
<ul>
<li>The locus passes QC filters (&quot;called&quot;)</li>
<li>The locus fails QC filters (&quot;filtered&quot;)</li>
<li>The locus either passes or fails QC filters (&quot;raw&quot;)</li>
</ul>
<h3>FunctionalClass</h3>
<p>FunctionalClass is a four-state stratification:</p>
<ul>
<li>The locus is a synonymous site (&quot;silent&quot;)</li>
<li>The locus is a missense site (&quot;missense&quot;)</li>
<li>The locus is a nonsense site (&quot;nonsense&quot;)</li>
<li>The locus is of any functional class (&quot;any&quot;)</li>
</ul>
<h3>CompRod</h3>
<p>CompRod is an N-state stratification, where N is the number of comp tracks bound to VariantEval.</p>
<h3>Degeneracy</h3>
<p>Degeneracy is a six-state stratification:</p>
<ul>
<li>The underlying base position in the codon is 1-fold degenerate (&quot;1-fold&quot;)</li>
<li>The underlying base position in the codon is 2-fold degenerate (&quot;2-fold&quot;)</li>
<li>The underlying base position in the codon is 3-fold degenerate (&quot;3-fold&quot;)</li>
<li>The underlying base position in the codon is 4-fold degenerate (&quot;4-fold&quot;)</li>
<li>The underlying base position in the codon is 6-fold degenerate (&quot;6-fold&quot;)</li>
<li>The underlying base position in the codon is degenerate at any level (&quot;all&quot;)</li>
</ul>
<p>See the [<a href="http://en.wikipedia.org/wiki/Genetic_code#Degeneracy">http://en.wikipedia.org/wiki/Genetic_code#Degeneracy</a> Wikipedia page on degeneracy] for more information.</p>
<h3>JexlExpression</h3>
<p>JexlExpression is an N-state stratification, where N is the number of JEXL expressions supplied to VariantEval.  See [[Using JEXL expressions]]</p>
<h3>Novelty</h3>
<p>Novelty is a three-state stratification:</p>
<ul>
<li>The locus overlaps the knowns comp track (usually the dbSNP track) (&quot;known&quot;)</li>
<li>The locus does not overlap the knowns comp track (&quot;novel&quot;)</li>
<li>The locus either overlaps or does not overlap the knowns comp track (&quot;all&quot;)</li>
</ul>
<h3>CountVariants</h3>
<p>CountVariants is an evaluation module that computes the following metrics:</p>
<table class="table table-striped">
<thead>
<tr>
<th style="text-align: left;">Metric</th>
<th style="text-align: left;">Definition</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">nProcessedLoci</td>
<td style="text-align: left;">Number of processed loci</td>
</tr>
<tr>
<td style="text-align: left;">nCalledLoci</td>
<td style="text-align: left;">Number of called loci</td>
</tr>
<tr>
<td style="text-align: left;">nRefLoci</td>
<td style="text-align: left;">Number of reference loci</td>
</tr>
<tr>
<td style="text-align: left;">nVariantLoci</td>
<td style="text-align: left;">Number of variant loci</td>
</tr>
<tr>
<td style="text-align: left;">variantRate</td>
<td style="text-align: left;">Variants per loci rate</td>
</tr>
<tr>
<td style="text-align: left;">variantRatePerBp</td>
<td style="text-align: left;">Number of variants per base</td>
</tr>
<tr>
<td style="text-align: left;">nSNPs</td>
<td style="text-align: left;">Number of snp loci</td>
</tr>
<tr>
<td style="text-align: left;">nInsertions</td>
<td style="text-align: left;">Number of insertion</td>
</tr>
<tr>
<td style="text-align: left;">nDeletions</td>
<td style="text-align: left;">Number of deletions</td>
</tr>
<tr>
<td style="text-align: left;">nComplex</td>
<td style="text-align: left;">Number of complex loci</td>
</tr>
<tr>
<td style="text-align: left;">nNoCalls</td>
<td style="text-align: left;">Number of no calls loci</td>
</tr>
<tr>
<td style="text-align: left;">nHets</td>
<td style="text-align: left;">Number of het loci</td>
</tr>
<tr>
<td style="text-align: left;">nHomRef</td>
<td style="text-align: left;">Number of hom ref loci</td>
</tr>
<tr>
<td style="text-align: left;">nHomVar</td>
<td style="text-align: left;">Number of hom var loci</td>
</tr>
<tr>
<td style="text-align: left;">nSingletons</td>
<td style="text-align: left;">Number of singletons</td>
</tr>
<tr>
<td style="text-align: left;">heterozygosity</td>
<td style="text-align: left;">heterozygosity per locus rate</td>
</tr>
<tr>
<td style="text-align: left;">heterozygosityPerBp</td>
<td style="text-align: left;">heterozygosity per base pair</td>
</tr>
<tr>
<td style="text-align: left;">hetHomRatio</td>
<td style="text-align: left;">heterozygosity to homozygosity ratio</td>
</tr>
<tr>
<td style="text-align: left;">indelRate</td>
<td style="text-align: left;">indel rate (insertion count + deletion count)</td>
</tr>
<tr>
<td style="text-align: left;">indelRatePerBp</td>
<td style="text-align: left;">indel rate per base pair</td>
</tr>
<tr>
<td style="text-align: left;">deletionInsertionRatio</td>
<td style="text-align: left;">deletion to insertion ratio</td>
</tr>
</tbody>
</table>
<h3>CompOverlap</h3>
<p>CompOverlap is an evaluation module that computes the following metrics:</p>
<table class="table table-striped">
<thead>
<tr>
<th style="text-align: left;">Metric</th>
<th style="text-align: left;">Definition</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">nEvalSNPs</td>
<td style="text-align: left;">number of eval SNP sites</td>
</tr>
<tr>
<td style="text-align: left;">nCompSNPs</td>
<td style="text-align: left;">number of comp SNP sites</td>
</tr>
<tr>
<td style="text-align: left;">novelSites</td>
<td style="text-align: left;">number of eval sites outside of comp sites</td>
</tr>
<tr>
<td style="text-align: left;">nVariantsAtComp</td>
<td style="text-align: left;">number of eval sites at comp sites (that is, sharing the same locus as a variant in the comp track, regardless of whether the alternate allele is the same)</td>
</tr>
<tr>
<td style="text-align: left;">compRate</td>
<td style="text-align: left;">percentage of eval sites at comp sites</td>
</tr>
<tr>
<td style="text-align: left;">nConcordant</td>
<td style="text-align: left;">number of concordant sites (that is, for the sites that share the same locus as a variant in the comp track, those that have the same alternate allele)</td>
</tr>
<tr>
<td style="text-align: left;">concordantRate</td>
<td style="text-align: left;">the concordance rate</td>
</tr>
</tbody>
</table>
<h4>Understanding the output of CompOverlap</h4>
<p>A SNP in the detection set is said to be 'concordant' if the position exactly matches an entry in dbSNP and the allele is the same.  To understand this and other output of CompOverlap, we shall examine a detailed example.  First, consider a fake dbSNP file (headers are suppressed so that one can see the important things):</p>
<pre><code class="pre_md"> $ grep -v '##' dbsnp.vcf
 #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
 1       10327   rs112750067     T       C       .       .       ASP;R5;VC=SNP;VP=050000020005000000000100;WGT=1;dbSNPBuildID=132</code class="pre_md"></pre>
<p>Now, a detection set file with a single sample, where the variant allele is the same as listed in dbSNP:</p>
<pre><code class="pre_md"> $ grep -v '##' eval_correct_allele.vcf
 #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT            001-6
 1       10327   .       T       C       5168.52 PASS    ...     GT:AD:DP:GQ:PL    0/1:357,238:373:99:3959,0,4059</code class="pre_md"></pre>
<p>Finally, a detection set file with a single sample, but the alternate allele differs from that in dbSNP:</p>
<pre><code class="pre_md"> $ grep -v '##' eval_incorrect_allele.vcf
 #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT            001-6
 1       10327   .       T       A       5168.52 PASS    ...     GT:AD:DP:GQ:PL    0/1:357,238:373:99:3959,0,4059</code class="pre_md"></pre>
<p>Running VariantEval with just the CompOverlap module:</p>
<pre><code class="pre_md"> $ java -jar $STING_DIR/dist/GenomeAnalysisTK.jar -T VariantEval \
        -R /seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta \
        -L 1:10327 \
        -B:dbsnp,VCF dbsnp.vcf \
        -B:eval_correct_allele,VCF eval_correct_allele.vcf \
        -B:eval_incorrect_allele,VCF eval_incorrect_allele.vcf \
        -noEV \
        -EV CompOverlap \
        -o eval.table</code class="pre_md"></pre>
<p>We find that the eval.table file contains the following:</p>
<pre><code class="pre_md"> $ grep -v '##' eval.table | column -t 
 CompOverlap  CompRod  EvalRod                JexlExpression  Novelty  nEvalVariants  nCompVariants  novelSites  nVariantsAtComp  compRate      nConcordant  concordantRate
 CompOverlap  dbsnp    eval_correct_allele    none            all      1              1              0           1                100.00000000  1            100.00000000
 CompOverlap  dbsnp    eval_correct_allele    none            known    1              1              0           1                100.00000000  1            100.00000000
 CompOverlap  dbsnp    eval_correct_allele    none            novel    0              0              0           0                0.00000000    0            0.00000000
 CompOverlap  dbsnp    eval_incorrect_allele  none            all      1              1              0           1                100.00000000  0            0.00000000
 CompOverlap  dbsnp    eval_incorrect_allele  none            known    1              1              0           1                100.00000000  0            0.00000000
 CompOverlap  dbsnp    eval_incorrect_allele  none            novel    0              0              0           0                0.00000000    0            0.00000000</code class="pre_md"></pre>
<p>As you can see, the detection set variant was listed under nVariantsAtComp (meaning the variant was seen at a position listed in dbSNP), but only the eval_correct_allele dataset is shown to be concordant at that site, because the allele listed in this dataset and dbSNP match.</p>
<h3>TiTvVariantEvaluator</h3>
<p>TiTvVariantEvaluator is an evaluation module that computes the following metrics:</p>
<table class="table table-striped">
<thead>
<tr>
<th style="text-align: left;">Metric</th>
<th style="text-align: left;">Definition</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">nTi</td>
<td style="text-align: left;">number of transition loci</td>
</tr>
<tr>
<td style="text-align: left;">nTv</td>
<td style="text-align: left;">number of transversion loci</td>
</tr>
<tr>
<td style="text-align: left;">tiTvRatio</td>
<td style="text-align: left;">the transition to transversion ratio</td>
</tr>
<tr>
<td style="text-align: left;">nTiInComp</td>
<td style="text-align: left;">number of comp transition sites</td>
</tr>
<tr>
<td style="text-align: left;">nTvInComp</td>
<td style="text-align: left;">number of comp transversion sites</td>
</tr>
<tr>
<td style="text-align: left;">TiTvRatioStandard</td>
<td style="text-align: left;">the transition to transversion ratio for comp sites</td>
</tr>
</tbody>
</table>