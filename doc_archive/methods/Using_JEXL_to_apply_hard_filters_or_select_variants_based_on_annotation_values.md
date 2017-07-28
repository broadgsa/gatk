## Using JEXL to apply hard filters or select variants based on annotation values

http://gatkforums.broadinstitute.org/gatk/discussion/1255/using-jexl-to-apply-hard-filters-or-select-variants-based-on-annotation-values

<h3>1. JEXL in a nutshell</h3>
<p>JEXL stands for Java EXpression Language. It's not a part of the GATK as such; it's a software library that can be used by Java-based programs like the GATK. It can be used for many things, but in the context of the GATK, it has one very specific use: making it possible to operate on subsets of variants from VCF files based on one or more annotations, using a single command. This is typically done with walkers such as <a href="https://www.broadinstitute.org/gatk/documentation/tooldocs/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php">VariantFiltration</a> and <a href="https://www.broadinstitute.org/gatk/documentation/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php">SelectVariants</a>.</p>
<hr />
<h3>2. Basic structure of JEXL expressions for use with the GATK</h3>
<p>In this context, a JEXL expression is a string (in the computing sense, <em>i.e.</em> a series of characters) that tells the GATK which annotations to look at and what selection rules to apply. </p>
<p>JEXL expressions contain three basic components: keys and values, connected by operators. For example, in this simple JEXL expression which selects variants whose quality score is greater than 30:</p>
<pre><code class="pre_md">"QUAL &gt; 30.0"</code class="pre_md"></pre>
<ul>
<li><code>QUAL</code> is a <strong>key</strong>: the name of the annotation we want to look at</li>
<li><code>30.0</code> is a <strong>value</strong>: the threshold that we want to use to evaluate variant quality against</li>
<li><code>&gt;</code> is an <strong>operator</strong>: it determines which &quot;side&quot; of the threshold we want to select</li>
</ul>
<p>The complete expression must be framed by double quotes. Within this, keys are strings (typically written in uppercase or CamelCase), and values can be either strings, numbers or booleans (TRUE or FALSE) -- but if they are strings the values must be framed by single quotes, as in the following example:</p>
<pre><code class="pre_md">"MY_STRING_KEY == 'foo'"</code class="pre_md"></pre>
<hr />
<h3>3. Evaluation on multiple annotations</h3>
<p>You can build expressions that calculate a metric based on two separate annotations, for example if you want to select variants for which quality (QUAL) divided by depth of coverage (DP) is below a certain threshold value:  </p>
<pre><code class="pre_md">"QUAL / DP &lt; 10.0"</code class="pre_md"></pre>
<p>You can also join multiple conditional statements with logical operators, for example if you want to select variants that have both sufficient quality (QUAL) <strong>and</strong> a certain depth of coverage (DP):</p>
<pre><code class="pre_md">"QUAL &gt; 30.0 &amp;&amp; DP == 10"</code class="pre_md"></pre>
<p>where <code>&amp;&amp;</code> is the logical &quot;AND&quot;.</p>
<p>Or if you want to select variants that have at least one of several conditions fulfilled:</p>
<pre><code class="pre_md">"QD &lt; 2.0 || ReadPosRankSum &lt; -20.0 || FS &gt; 200.0"</code class="pre_md"></pre>
<p>where <code>||</code> is the logical &quot;OR&quot;.</p>
<hr />
<h3>4. Filtering on sample/genotype-level properties</h3>
<p>You can also filter individual samples/genotypes in a VCF based on information from the FORMAT field. Variant Filtration will add the sample-level FT tag to the FORMAT field of filtered samples. Note however that this does not affect the record's FILTER tag. This is still a work in progress and isn't quite as flexible and powerful yet as we'd like it to be. For now, you can filter based on most fields as normal (e.g. GQ &lt; 5.0), but the GT (genotype) field is an exception. We have put in convenience methods to enable filtering out heterozygous calls (isHet == 1), homozygous-reference calls (isHomRef == 1), and homozygous-variant calls (isHomVar == 1).</p>
<hr />
<h3>5. Important caveats</h3>
<h4>Sensitivity to case and type</h4>
<p>You're probably used to case being important (whether letters are lowercase or UPPERCASE) but now you need to also pay attention to the type of value that is involved -- for example, numbers are differentiated between integers and floats (essentially, non-integers). These points are especially important to keep in mind:</p>
<ul>
<li>Case</li>
</ul>
<p>Currently, VCF INFO field keys are case-sensitive. That means that if you have a <code>QUAL</code> field in uppercase in your VCF record, the system will not recognize it if you write it differently (<code>Qual</code>, <code>qual</code> or whatever) in your JEXL expression. </p>
<ul>
<li>Type</li>
</ul>
<p>The types (<em>i.e.</em> string, integer, non-integer or boolean) used in your expression must be exactly the same as that of the value you are trying to evaluate. In other words, if you have a QUAL field with non-integer values (<em>e.g.</em> <strong>45.3</strong>) and your filter expression is written as an integer (<em>e.g.</em> &quot;QUAL &lt; <strong>50</strong>&quot;), the system will throw a hissy fit (<em>aka</em> a Java exception).</p>
<h4>Complex queries</h4>
<p>We highly recommend that complex expressions involving multiple AND/OR operations be split up into separate expressions whenever  possible to avoid confusion. If you are using complex expressions, make sure to test them on a panel of different sites with several combinations of yes/no criteria.</p>
<hr />
<h3>6. More complex JEXL magic</h3>
<p>Note that this last part is fairly advanced and not for the faint of heart. To be frank, it's also explained rather more briefly than the topic deserves. But if there's enough demand for this level of usage (click the &quot;view in forum&quot; link and leave a comment) we'll consider producing a full-length tutorial. </p>
<h4>Introducing the VariantContext object</h4>
<p>When you use SelectVariants with JEXL, what happens under the hood is that the program accesses something called the VariantContext, which is a representation of the variant call with all its annotation information. The VariantContext is technically not part of GATK; it's part of the <code>variant</code> library included within the Picard tools source code, which GATK uses for convenience. </p>
<p>The reason we're telling you about this is that you can actually make more complex queries than what the GATK offers convenience functions for, provided you're willing to do a little digging into the VariantContext methods. This will allow you to leverage the full range of capabilities of the underlying objects from the command line.</p>
<p>In a nutshell, the VariantContext is available through the <code>vc</code> variable, and you just need to add method calls to that variable in your command line. The bets way to find out what methods are available is to read the <a href="http://sourceforge.net/p/picard/code/HEAD/tree/trunk/src/java/org/broadinstitute/variant/variantcontext/VariantContext.java">VariantContext documentation</a> on the Picard tools source code repository (on SourceForge), but we list a few examples below to whet your appetite.</p>
<h4>Using VariantContext directly</h4>
<p>For example, suppose I want to use SelectVariants to select all of the sites where sample NA12878 is homozygous-reference. This can be accomplished by assessing the underlying VariantContext as follows:</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T SelectVariants -R b37/human_g1k_v37.fasta --variant my.vcf -select 'vc.getGenotype("NA12878").isHomRef()'</code class="pre_md"></pre>
<p>Groovy, right? Now here's a more sophisticated example of JEXL expression that finds all novel variants in the total set with allele frequency &gt; 0.25 but not 1, is not filtered, and is non-reference in 01-0263 sample:</p>
<pre><code class="pre_md">! vc.getGenotype("01-0263").isHomRef() &amp;&amp; (vc.getID() == null || vc.getID().equals(".")) &amp;&amp; AF &gt; 0.25 &amp;&amp; AF &lt; 1.0 &amp;&amp; vc.isNotFiltered() &amp;&amp; vc.isSNP() -o 01-0263.high_freq_novels.vcf -sn 01-0263</code class="pre_md"></pre>
<h4>Using the VariantContext to evaluate boolean values</h4>
<p>The classic way of evaluating a boolean goes like this:</p>
<pre><code class="pre_md">java -Xmx4g -jar GenomeAnalysisTK.jar -T SelectVariants -R b37/human_g1k_v37.fasta --variant my.vcf -select 'DB'</code class="pre_md"></pre>
<p>But you can also use the VariantContext object like this:</p>
<pre><code class="pre_md">java -Xmx4g -jar GenomeAnalysisTK.jar -T SelectVariants -R b37/human_g1k_v37.fasta --variant my.vcf -select 'vc.hasAttribute("DB")'</code class="pre_md"></pre>
<h4>Using VariantContext to access annotations in multiallelic sites</h4>
<p>The order of alleles in the VariantContext object is not guaranteed to be the same as in the VCF output, so accessing the AF by an index derived from a scrambled alleles array is dangerous. However! If we have the sample genotypes, there's a workaround:</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T SelectVariants -R reference.fasta -V multiallelics.vcf -select 'vc.hasGenotypes() &amp;&amp; vc.getCalledChrCount(vc.getAltAlleleWithHighestAlleleCount())/(1.0*vc.getCalledChrCount()) &gt; 0.1' -o multiHighAC.vcf</code class="pre_md"></pre>
<p>The odd 1.0 is there because otherwise we're dividing two integers, which will always yield 0. The <code>vc.hasGenotypes()</code> is extra error checking. This might be slow for large files, but we could use something like this if performance is a concern:</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T SelectVariants -R reference.fasta -V multiallelics.vcf -select 'vc.isBiallelic() ? AF &gt; 0.1 : vc.hasGenotypes() &amp;&amp; vc.getCalledChrCount(vc.getAltAlleleWithHighestAlleleCount())/(1.0*vc.getCalledChrCount()) &gt; 0.1' -o multiHighAC.vcf</code class="pre_md"></pre>
<p>Where hopefully the ternary expression shortcuts the extra <code>vc</code> calls for all the biallelics.</p>
<h4>Using JEXL to evaluate arrays</h4>
<p>Sometimes you might want to write a JEXL expression to evaluate e.g. the AD (allelic depth) field in the FORMAT column.  However, the AD is technically not an integer; rather it is a list (array) of integers.  One can evaluate the array data using the &quot;.&quot; operator.   Here's an example:</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T SelectVariants -R b37/human_g1k_v37.fasta --variant my.vcf -select 'vc.getGenotype("NA12878").getAD().0 &gt; 10'</code class="pre_md"></pre>
<p>If you would like to select sites where the alternate allele frequency is greater than 50%, you can use the following expression:</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T SelectVariants -R b37/human_g1k_v37.fasta --variant my.vcf -select vc.getGenotype("NA12878").getAD().1 / vc.getGenotype("NA12878").getDP() &gt; 0.50</code class="pre_md"></pre>