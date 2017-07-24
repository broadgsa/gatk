## (howto) Run the genotype refinement workflow

http://gatkforums.broadinstitute.org/gatk/discussion/4727/howto-run-the-genotype-refinement-workflow

<h3>Overview</h3>
<p>This tutorial describes step-by-step instruction for applying the Genotype Refinement workflow (described in <a href="https://www.broadinstitute.org/gatk/guide/article?id=4723">this method article</a>) to your data.</p>
<hr />
<h3>Step 1: Derive posterior probabilities of genotypes</h3>
<p>In this first step, we are deriving the posteriors of genotype calls in our callset, <code>recalibratedVariants.vcf</code>, which just came out of the VQSR filtering step; it contains among other samples a trio of individuals (mother, father and child) whose family structure is described in the pedigree file <code>trio.ped</code> (which you need to supply). To do this, we are using the most comprehensive set of high confidence SNPs available to us, a set of sites from Phase 3 of the 1000 Genomes project (available in our resource bundle), which we pass via the <code>--supporting</code> argument.</p>
<pre><code class="pre_md"> java -jar GenomeAnalysisToolkit.jar -R human_g1k_v37_decoy.fasta -T CalculateGenotypePosteriors --supporting 1000G_phase3_v4_20130502.sites.vcf -ped trio.ped -V recalibratedVariants.vcf -o recalibratedVariants.postCGP.vcf</code class="pre_md"></pre>
<p>This produces the output file <code>recalibratedVariants.postCGP.vcf</code>, in which the posteriors have been annotated wherever possible.</p>
<hr />
<h3>Step 2: Filter low quality genotypes</h3>
<p>In this second, very simple step, we are tagging low quality genotypes so we know not to use them in our downstream analyses. We use Q20 as threshold for quality, which means that any passing genotype has a 99% chance of being correct. </p>
<pre><code class="pre_md">java -jar $GATKjar -T VariantFiltration -R $bundlePath/b37/human_g1k_v37_decoy.fasta -V recalibratedVariants.postCGP.vcf -G_filter "GQ &lt; 20.0" -G_filterName lowGQ -o recalibratedVariants.postCGP.Gfiltered.vcf</code class="pre_md"></pre>
<p>Note that in the resulting VCF, the genotypes that failed the filter are still present, but they are tagged <code>lowGQ</code> with the FT tag of the FORMAT field.</p>
<hr />
<h3>Step 3: Annotate possible <em>de novo</em> mutations</h3>
<p>In this third and final step, we tag variants for which at least one family in the callset shows evidence of a <em>de novo</em> mutation based on the genotypes of the family members. </p>
<pre><code class="pre_md">java -jar $GATKjar -T VariantAnnotator -R $bundlePath/b37/human_g1k_v37_decoy.fasta -V recalibratedVariants.postCGP.Gfiltered.vcf -A PossibleDeNovo -ped trio.ped -o recalibratedVariants.postCGP.Gfiltered.deNovos.vcf</code class="pre_md"></pre>
<p>The annotation output will include a list of the children with possible <em>de novo</em> mutations, classified as either high or low confidence.</p>
<p>See section 3 of the method article for a complete description of annotation outputs and section 4 for an example of a call and the interpretation of the annotation values.</p>