## Writing and working with reference metadata classes

http://gatkforums.broadinstitute.org/gatk/discussion/1350/writing-and-working-with-reference-metadata-classes

<h2>Brief introduction to reference metadata (RMDs)</h2>
<p><em>Note that the <code>-B</code> flag referred to below is deprecated; these docs need to be updated</em></p>
<p>The GATK allows you to process arbitrary numbers of reference metadata (RMD) files inside of walkers (previously we called this reference ordered data, or ROD).  Common RMDs are things like dbSNP, VCF call files, and refseq annotations.  The only real constraints on RMD files is that:</p>
<ul>
<li>
<p>They must contain information necessary to provide contig and position data for each element to the GATK engine so it knows with what loci to associate the RMD element.</p>
</li>
<li>
<p>The file must be sorted with regard to the reference fasta file so that data can be accessed sequentially by the engine.</p>
</li>
<li>The file must have a <a href="http://gatkforums.broadinstitute.org/discussion/1349/tribble">Tribble</a> RMD parsing class associated with the file type so that elements in the RMD file can be parsed by the engine.</li>
</ul>
<p>Inside of the GATK the RMD system has the concept of RMD tracks, which associate an arbitrary string name with the data in the associated RMD file.  For example, the <code>VariantEval</code> module uses the named track <code>eval</code> to get calls for evaluation, and <code>dbsnp</code> as the track containing the database of known variants.</p>
<h2>How do I get reference metadata files into my walker?</h2>
<p>RMD files are extremely easy to get into the GATK using the <code>-B</code> syntax:</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -R Homo_sapiens_assembly18.fasta -T PrintRODs -B:variant,VCF calls.vcf</code class="pre_md"></pre>
<p>In this example, the GATK will attempt to parse the file <code>calls.vcf</code> using the VCF parser and bind the VCF data to the RMD track named <code>variant</code>.</p>
<p>In general, you can provide as many RMD bindings to the GATK as you like:</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -R Homo_sapiens_assembly18.fasta -T PrintRODs -B:calls1,VCF calls1.vcf -B:calls2,VCF calls2.vcf</code class="pre_md"></pre>
<p>Works just as well.  Some modules may require specifically named RMD tracks -- like <code>eval</code> above -- and some are happy to just assess all RMD tracks of a certain class and work with those -- like <code>VariantsToVCF</code>.</p>
<h3>1. Directly getting access to a single named track</h3>
<p>In this snippet from <code>SNPDensityWalker</code>, we grab the <code>eval</code> track as a <code>VariantContext</code> object, only for the variants that are of type SNP:</p>
<pre><code class="pre_md">public Pair&lt;VariantContext, GenomeLoc&gt; map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
    VariantContext vc = tracker.getVariantContext(ref, "eval", EnumSet.of(VariantContext.Type.SNP), context.getLocation(), false);
}</code class="pre_md"></pre>
<h3>2. Grabbing anything that's convertable to a VariantContext</h3>
<p>From <code>VariantsToVCF</code> we call the helper function <code>tracker.getVariantContexts</code> to look at all of the RMDs and convert what it can to <code>VariantContext</code> objects.</p>
<pre><code class="pre_md">Allele refAllele = new Allele(Character.toString(ref.getBase()), true);
Collection&lt;VariantContext&gt; contexts = tracker.getVariantContexts(INPUT_RMD_NAME, ALLOWED_VARIANT_CONTEXT_TYPES, context.getLocation(), refAllele, true, false);</code class="pre_md"></pre>
<h3>3. Looking at all of the RMDs</h3>
<p>Here's a totally general code snippet from <code>PileupWalker.java</code>.  This code, as you can see, iterates over all of the GATKFeature objects in the reference ordered data, converting each RMD to a string and capturing these strings in a list.  It finally grabs the dbSNP binding specifically for a more detailed string conversion, and then binds them all up in a single string for display along with the read pileup.</p>
<p>private String getReferenceOrderedData( RefMetaDataTracker tracker ) {
ArrayList<String> rodStrings = new ArrayList<String>();
for ( GATKFeature datum : tracker.getAllRods() ) {
if ( datum != null &amp;&amp; ! (datum.getUnderlyingObject() instanceof DbSNPFeature)) {
rodStrings.add(((ReferenceOrderedDatum)datum.getUnderlyingObject()).toSimpleString()); // TODO: Aaron: this line still survives, try to remove it
}
}
String rodString = Utils.join(&quot;, &quot;, rodStrings);</p>
<pre><code class="pre_md">        DbSNPFeature dbsnp = tracker.lookup(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME, DbSNPFeature.class);

        if ( dbsnp != null)
            rodString += DbSNPHelper.toMediumString(dbsnp);

        if ( !rodString.equals("") )
            rodString = "[ROD: " + rodString + "]";

        return rodString;
}</code class="pre_md"></pre>
<h2>How do I write my own RMD types?</h2>
<p>Tracks of reference metadata are loaded using the <a href="http://gatkforums.broadinstitute.org/discussion/1349/tribble">Tribble</a> infrastructure.  Tracks are loaded using the feature codec and underlying type information.  See the <a href="http://gatkforums.broadinstitute.org/discussion/1349/tribble">Tribble documentation</a> for more information.</p>
<p>Tribble codecs that are in the classpath are automatically found; the GATK discovers all classes that implement the <code>FeatureCodec</code> class. Name resolution occurs using the <code>-B</code> type parameter, i.e. if the user specified:   </p>
<pre><code class="pre_md">-B:calls1,VCF calls1.vcf</code class="pre_md"></pre>
<p>The GATK looks for a <code>FeatureCodec</code> called <code>VCFCodec.java</code> to decode the record type.  Alternately, if the user specified:</p>
<pre><code class="pre_md">-B:calls1,MYAwesomeFormat calls1.maft</code class="pre_md"></pre>
<p>THe GATK would look for a codec called <code>MYAwesomeFormatCodec.java</code>.  This look-up is not case sensitive, i.e. it will resolve <code>MyAwEsOmEfOrMaT</code> as well, though why you would want to write something so painfully ugly to read is beyond us.</p>