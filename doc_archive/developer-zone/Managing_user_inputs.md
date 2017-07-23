## Managing user inputs

http://gatkforums.broadinstitute.org/gatk/discussion/1325/managing-user-inputs

<h3>1. Naming walkers</h3>
<p>Users identify which GATK walker to run by specifying a walker name via the <code>--analysis_type</code> command-line argument.  By default, the GATK will derive the walker name from a walker by taking the name of the walker class and removing packaging information from the start of the name, and removing the trailing text <code>Walker</code> from the end of the name, if it exists.  For example, the GATK would, by default, assign the name <code>PrintReads</code> to the walker class <code>org.broadinstitute.sting.gatk.walkers.PrintReadsWalker</code>.  To override the default walker name, annotate your walker class with <code>@WalkerName("&lt;my name&gt;")</code>.</p>
<h3>2. Requiring / allowing primary inputs</h3>
<p>Walkers can flag exactly which primary data sources are allowed and required for a given walker.   Reads, the reference, and reference-ordered data are currently considered primary data sources.  Different traversal types have different default requirements for reads and reference, but currently no traversal types require reference-ordered data by default.  You can add requirements to your walker with the <code>@Requires</code> / <code>@Allows</code> annotations as follows:</p>
<pre><code class="pre_md">@Requires(DataSource.READS)
@Requires({DataSource.READS,DataSource.REFERENCE})
@Requires(value={DataSource.READS,DataSource.REFERENCE})
@Requires(value=DataSource.REFERENCE})</code class="pre_md"></pre>
<p>By default, all parameters are allowed unless you lock them down with the <code>@Allows</code> attribute.  The command:</p>
<pre><code class="pre_md">@Allows(value={DataSource.READS,DataSource.REFERENCE})</code class="pre_md"></pre>
<p>will only allow the reads and the reference. Any other primary data sources will cause the system to exit with an error.  </p>
<p>Note that as of August 2011, the GATK no longer supports RMD the <code>@Requires</code> and <code>@Allows</code> syntax, as these have moved to the standard <code>@Argument</code> system.</p>
<h3>3. Command-line argument tagging</h3>
<p>Any command-line argument can be tagged with a comma-separated list of freeform tags.  </p>
<p>The syntax for tags is as follows:</p>
<pre><code class="pre_md">-&lt;argument&gt;:&lt;tag1&gt;,&lt;tag2&gt;,&lt;tag3&gt; &lt;argument value&gt;</code class="pre_md"></pre>
<p>for example:</p>
<pre><code class="pre_md">-I:tumor &lt;my tumor data&gt;.bam
-eval,VCF yri.trio.chr1.vcf</code class="pre_md"></pre>
<p>There is currently no mechanism in the GATK to validate either the number of tags supplied or the content of those tags.</p>
<p>Tags can be accessed from within a walker by calling <code>getToolkit().getTags(argumentValue)</code>, where <code>argumentValue</code> is the
parsed contents of the command-line argument to inspect.</p>
<h4>Applications</h4>
<p>The GATK currently has comprehensive support for tags on two built-in argument types:</p>
<ul>
<li>
<p><code>-I,--input_file &lt;input_file&gt;</code></p>
<p>Input BAM files and BAM file lists can be tagged with any type.  When a BAM file list is tagged, the tag is applied to each listed BAM file.  </p>
</li>
</ul>
<p>From within a walker, use the following code to access the supplied tag or tags:</p>
<pre><code class="pre_md">getToolkit().getReaderIDForRead(read).getTags();</code class="pre_md"></pre>
<ul>
<li>
<p>Input RODs, e.g. `-V <rod>' or '-eval <rod>'</p>
<p>Tags are used to specify ROD name and ROD type.  There is currently no support for adding additional tags.  See the ROD system documentation for more details.</p>
</li>
</ul>
<h3>4. Adding additional command-line arguments</h3>
<p>Users can create command-line arguments for walkers by creating public member variables annotated with <code>@Argument</code> in the walker. The <code>@Argument</code> annotation takes a number of differentparameters:</p>
<ul>
<li>
<p><code>fullName</code></p>
<p>The full name of this argument. Defaults to the <code>toLowerCase()</code>’d member name. When specifying <code>fullName</code> on the command line, prefix with a double dash (<code>--</code>).</p>
</li>
<li>
<p><code>shortName</code> </p>
<p>The alternate, short name for this argument. Defaults to the first letter of the member name.  When specifying shortName on the command line, prefix with a single dash (<code>-</code>).</p>
</li>
<li>
<p><code>doc</code> </p>
<p>Documentation for this argument. Will appear in help output when a user either requests help with the –-help (-h) argument or when a user specifies an invalid set of arguments.  Documentation is the only argument that is always required.</p>
</li>
<li>
<p><code>required</code> </p>
<p>Whether the argument is required when used with this walker. Default is <code>required = true</code>.</p>
</li>
<li>
<p><code>exclusiveOf</code> </p>
<p>Specifies that this argument is mutually exclusive of another argument in the same walker.  Defaults to not mutually exclusive of any other arguments.</p>
</li>
<li>
<p><code>validation</code> </p>
<p>Specifies a regular expression used to validate the contents of the command-line argument.  If the text provided by the user does not match this regex, the GATK will abort with an error.</p>
</li>
</ul>
<p>By default, all command-line arguments will appear in the help system.  To prevent new and debugging arguments from appearing in the help system,
you can add the <code>@Hidden</code> tag below the <code>@Argument</code> annotation, hiding it from the help system but allowing users to supply it on the command-line.
Please use this functionality sparingly to avoid walkers with hidden command-line options that are required for production use.</p>
<h4>Passing Command-Line Arguments</h4>
<p>Arguments can be passed to the walker using either the full name or the short name. If passing arguments using the full name, the syntax is <code>−−&lt;arg full name&gt; &lt;value&gt;</code>.</p>
<pre><code class="pre_md">--myint 6</code class="pre_md"></pre>
<p>If passing arguments using the short name, the syntax is <code>-&lt;arg short name&gt; &lt;value&gt;</code>. Note that there is a space between the short name and the value:</p>
<pre><code class="pre_md">-m 6</code class="pre_md"></pre>
<p>Boolean (class) and boolean (primitive) arguments are a special in that they require no argument. The presence of a boolean indicates true, and its absence indicates false. The following example sets a flag to true.</p>
<pre><code class="pre_md">-B</code class="pre_md"></pre>
<h4>Supplemental command-line argument annotations</h4>
<p>Two additional annotations can influence the behavior of command-line arguments.</p>
<ul>
<li>
<p><code>@Hidden</code> </p>
<p>Adding this annotation to an @Argument tells the help system to avoid displaying any evidence that this argument exists.  This can be used to add additional debugging arguments that aren't suitable for mass consumption.</p>
</li>
<li>
<p><code>@Deprecated</code> </p>
<p>Forces the GATK to throw an exception if this argument is supplied on the command-line.  This can be used to supply extra documentation to the user as command-line parameters change for walkers that are in flux.</p>
</li>
</ul>
<h4>Examples</h4>
<p>Create an required int parameter with full name <code>–myint</code>, short name <code>-m</code>. Pass this argument by adding <code>–myint 6</code> or <code>-m 6</code> to the command line.</p>
<pre><code class="pre_md">import org.broadinstitute.sting.utils.cmdLine.Argument;
public class HelloWalker extends ReadWalker&lt;Integer,Long&gt; {
    @Argument(doc="my integer")
    public int myInt;</code class="pre_md"></pre>
<p>Create an optional float parameter with full name <code>–myFloatingPointArgument</code>, short name <code>-m</code>. Pass this argument by adding <code>–myFloatingPointArgument 2.71</code> or <code>-m 2.71</code>.</p>
<pre><code class="pre_md">import org.broadinstitute.sting.utils.cmdLine.Argument;
public class HelloWalker extends ReadWalker&lt;Integer,Long&gt; {
    @Argument(fullName="myFloatingPointArgument",doc="a floating point argument",required=false)
    public float myFloat;</code class="pre_md"></pre>
<p>The GATK will parse the argument differently depending on the type of the public member variable’s type. Many different argument types are supported, including primitives and their wrappers, arrays, typed and untyped collections, and any type with a String constructor. When the GATK cannot completely infer the type (such as in the case of untyped collections), it will assume that the argument is a String. GATK is aware of concrete implementations of some interfaces and abstract classes. If the argument’s member variable is of type <code>List</code> or <code>Set</code>, the GATK will fill the member variable with a concrete <code>ArrayList</code> or <code>TreeSet</code>, respectively. Maps are not currently supported.</p>
<h3>5. Additional argument types: @Input, @Output</h3>
<p>Besides <code>@Argument</code>, the GATK provides two additional types for command-line arguments: <code>@Input</code> and <code>@Output</code>.  These two inputs are very similar to <code>@Argument</code> but act as flags to indicate dataflow to <a href="http://gatkforums.broadinstitute.org/discussion/1306/overview-of-queue">Queue</a>, our pipeline management software.</p>
<ul>
<li>
<p>The <code>@Input</code> tag indicates that the contents of the tagged field represents a file that will be read by the walker.</p>
</li>
<li>The <code>@Output</code> tag indicates that the contents of the tagged field represents a file that will be written by the walker, for consumption by downstream walkers.</li>
</ul>
<p>We're still determining the best way to model walker dependencies in our pipeline.  As we determine best practices, we'll post them here.</p>
<h3>6. Getting access to Reference Ordered Data (RMD) with @Input and RodBinding<T></h3>
<p>As of August 2011, the GATK now provides a clean mechanism for creating walker <code>@Input</code> arguments and using these arguments to access <code>Reference Meta Data</code> provided by the <code>RefMetaDataTracker</code> in the <code>map()</code> call.  This mechanism is preferred to the old implicit string-based mechanism, which has been retired.</p>
<p>At a very high level, the new <code>RodBindings</code> provide a handle for a walker to obtain the <code>Feature</code> records from <code>Tribble</code> from a <code>map()</code> call, specific to a command line binding provided by the user.  This can be as simple as a single ROD file argument|one-to-one binding between a command line argument and a track, or as complex as an argument argument accepting multiple command line arguments, each with a specific name.  The <code>RodBindings</code> are generic and type specific, so you can require users to provide files that emit <code>VariantContext</code>s, <code>BedTable</code>s, etc, or simply the root type <code>Feature</code> from <code>Tribble</code>.   Critically, the <code>RodBindings</code> interact nicely with the GATKDocs system, so you can provide summary and detailed documentation for each <code>RodBinding</code> accepted by your walker.  </p>
<h4>A single ROD file argument</h4>
<p>Suppose you have a walker that uses a single track of <code>VariantContext</code>s, such as <code>SelectVariants</code>, in its calculation.  You declare a standard GATK-style <code>@Input</code> argument in the walker, of type <code>RodBinding&lt;VariantContext&gt;</code>: </p>
<pre><code class="pre_md">@Input(fullName="variant", shortName = "V", doc="Select variants from this VCF file", required=true)
public RodBinding&lt;VariantContext&gt; variants;</code class="pre_md"></pre>
<p>This will require the user to provide a command line option <code>--variant:vcf my.vcf</code> to your walker.  To get access to your variants, in the <code>map()</code> function you provide the variants variable to the tracker, as in:</p>
<pre><code class="pre_md">Collection&lt;VariantContext&gt; vcs = tracker.getValues(variants, context.getLocation());</code class="pre_md"></pre>
<p>which returns all of the <code>VariantContexts</code> in variants that start at <code>context.getLocation()</code>.  See <code>RefMetaDataTracker</code> in the javadocs to see the full range of getter routines.</p>
<p>Note that, as with regular tribble tracks, you have to provide the <code>Tribble</code> type of the file as a tag to the argument (<code>:vcf</code>).  The system now checks up front that the corresponding <code>Tribble</code> codec produces <code>Features</code> that are type-compatible with the type of the <code>RodBinding&lt;T&gt;</code>. </p>
<h4>RodBindings are generic</h4>
<p>The <code>RodBinding</code> class is generic, parameterized as <code>RodBinding&lt;T extends Feature&gt;</code>.  This <code>T</code> class describes the type of the <code>Feature</code> required by the walker.  The best practice for declaring a <code>RodBinding</code> is to choose the most general <code>Feature</code> type that will allow your walker to work.  For example, if all you really care about is whether a <code>Feature</code> overlaps the site in map, you can use <code>Feature</code> itself, which supports this, and will allow any <code>Tribble</code> type to be provided, using a <code>RodBinding&lt;Feature&gt;</code>.  If you are manipulating <code>VariantContext</code>s, you should declare a <code>RodBinding&lt;VariantContext&gt;</code>, which will restrict automatically the user to providing <code>Tribble</code> types that can create a object consistent with the <code>VariantContext</code> class (a <code>VariantContext</code> itself or subclass).</p>
<p>Note that in multi-argument <code>RodBindings</code>, as <code>List&lt;RodBinding&lt;T&gt;&gt;</code> arg, the system will require all files provided here to provide an object of type <code>T</code>.  So <code>List&lt;RodBinding&lt;VariantContext&gt;&gt;</code> arg requires all <code>-arg</code> command line arguments to bind to files that produce <code>VariantContext</code>s.</p>
<h4>An argument that can be provided any number of times</h4>
<p>The <code>RodBinding</code> system supports the standard <code>@Argument</code> style of allowing a <code>vararg</code> argument by wrapping it in a Java collection.  For example, if you want to allow users to provide any number of comp tracks to your walker, simply declare a <code>List&lt;RodBinding&lt;VariantContext&gt;&gt;</code> field:</p>
<pre><code class="pre_md">@Input(fullName="comp", shortName = "comp", doc="Comparison variants from this VCF file", required=true)
public List&lt;RodBinding&lt;VariantContext&gt;&gt; comps;</code class="pre_md"></pre>
<p>With this declaration, your walker will accept any number of <code>-comp</code> arguments, as in:</p>
<pre><code class="pre_md">-comp:vcf 1.vcf -comp:vcf 2.vcf -comp:vcf 3.vcf</code class="pre_md"></pre>
<p>For such a command line, the comps field would be initialized to the List with three <code>RodBindings</code>, the first binding to <code>1.vcf</code>, the second to <code>2.vcf</code> and finally the third to <code>3.vcf</code>.  </p>
<p>Because this is a required argument, at least one <code>-comp</code> must be provided.  <code>Vararg</code> <code>@Input</code> <code>RodBindings</code> can be optional, but you should follow proper <code>vararg</code>s style to get the best results.</p>
<h4>Proper handling of optional arguments</h4>
<p>If you want to make a RodBinding optional, you first need to tell the <code>@Input</code> argument that its options (<code>required=false</code>):</p>
<pre><code class="pre_md">@Input(fullName="discordance", required=false)
private RodBinding&lt;VariantContext&gt; discordanceTrack;</code class="pre_md"></pre>
<p>The GATK automagically sets this field to the value of the special static constructor method <code>makeUnbound(Class c)</code> to create a special &quot;unbound&quot; <code>RodBinding</code> here.  This unbound object is type safe, can be safely passed to the <code>RefMetaDataTracker</code> get methods, and is guaranteed to never return any values.  It also returns <code>false</code> when the <code>isBound()</code> method is called.</p>
<p>An example usage of <code>isBound</code> is to conditionally add header lines, as in:</p>
<pre><code class="pre_md">if ( mask.isBound() ) {
    hInfo.add(new VCFFilterHeaderLine(MASK_NAME, "Overlaps a user-input mask"));
}</code class="pre_md"></pre>
<p>The case for <code>vararg</code> style <code>RodBindings</code> is slightly different.  If you want, as above, users to be able to omit the <code>-comp</code> track entirely, you should initialize the value of the collection to the appropriate <code>emptyList</code>/<code>emptySet</code> in <code>Collections</code>:</p>
<pre><code class="pre_md">@Input(fullName="comp", shortName = "comp", doc="Comparison variants from this VCF file", required=false)
public List&lt;RodBinding&lt;VariantContext&gt;&gt; comps = Collections.emptyList();</code class="pre_md"></pre>
<p>which will ensure that <code>comps.isEmpty()</code> is true when no <code>-comp</code> is provided.</p>
<h4>Implicit and explicit names for RodBindings</h4>
<pre><code class="pre_md">@Input(fullName="variant", shortName = "V", doc="Select variants from this VCF file", required=true)
public RodBinding&lt;VariantContext&gt; variants;</code class="pre_md"></pre>
<p>By default, the <code>getName()</code> method in <code>RodBinding</code> returns the <code>fullName</code> of the <code>@Input</code>.  This can be overloaded on the command-line by providing not one but two tags.  The first tag is interpreted as the name for the binding, and the second as the type.  As in:</p>
<pre><code class="pre_md">-variant:vcf foo.vcf     =&gt; getName() == "variant"
-variant:foo,vcf foo.vcf =&gt; getName() == "foo"</code class="pre_md"></pre>
<p>This capability is useful when users need to provide more meaningful names for arguments, especially with variable arguments.  For example, in <code>VariantEval</code>, there's a <code>List&lt;RodBinding&lt;VariantContext&gt;&gt;</code> comps, which may be <code>dbsnp</code>, <code>hapmap</code>, etc.  This would be declared as:</p>
<pre><code class="pre_md">@Input(fullName="comp", shortName = "comp", doc="Comparison variants from this VCF file", required=true)
public List&lt;RodBinding&lt;VariantContext&gt;&gt; comps;</code class="pre_md"></pre>
<p>where a normal command line usage would look like:</p>
<pre><code class="pre_md">-comp:hapmap,vcf hapmap.vcf -comp:omni,vcf omni.vcf -comp:1000g,vcf 1000g.vcf</code class="pre_md"></pre>
<p>In the code, you might have a loop that looks like:</p>
<pre><code class="pre_md">for ( final RodBinding comp : comps )
    for ( final VariantContext vc : tracker.getValues(comp, context.getLocation())
        out.printf("%s has a binding at %s%n", comp.getName(), getToolkit().getGenomeLocParser.createGenomeLoc(vc)); </code class="pre_md"></pre>
<p>which would print out lines that included things like:</p>
<pre><code class="pre_md">hapmap has a binding at 1:10
omni has a binding at 1:20
hapmap has a binding at 1:30
1000g has a binding at 1:30</code class="pre_md"></pre>
<p>This last example begs the question -- what happens with <code>getName()</code> when explicit names are not provided?  The system goes out of its way to provide reasonable names for the variables: </p>
<ul>
<li>
<p>The first occurrence is named for the <code>fullName</code>, where <code>comp</code></p>
</li>
<li>Subsequent occurrences are postfixed with an integer count, starting at 2, so <code>comp2</code>, <code>comp3</code>, etc.</li>
</ul>
<p>In the above example, the command line </p>
<pre><code class="pre_md">-comp:vcf hapmap.vcf -comp:vcf omni.vcf -comp:vcf 1000g.vcf</code class="pre_md"></pre>
<p>would emit</p>
<pre><code class="pre_md">comp has a binding at 1:10
comp2 has a binding at 1:20
comp has a binding at 1:30
comp3 has a binding at 1:30</code class="pre_md"></pre>
<h4>Dynamic type resolution</h4>
<p>The new <code>RodBinding</code> system supports a simple form of dynamic type resolution.  If the input filetype can be specially associated with a single <code>Tribble</code> type (as VCF can), then you can omit the type entirely from the the command-line binding of a <code>RodBinding</code>!</p>
<p>So whereas a full command line would look like:</p>
<pre><code class="pre_md">-comp:hapmap,vcf hapmap.vcf -comp:omni,vcf omni.vcf -comp:1000g,vcf 1000g.vcf</code class="pre_md"></pre>
<p>because these are VCF files they could technically be provided as:</p>
<pre><code class="pre_md">-comp:hapmap hapmap.vcf -comp:omni omni.vcf -comp:1000g 1000g.vcf</code class="pre_md"></pre>
<p>If you don't care about naming, you can now say:</p>
<pre><code class="pre_md">-comp hapmap.vcf -comp omni.vcf -comp 1000g.vcf</code class="pre_md"></pre>
<h4>Best practice for documenting a RodBinding</h4>
<p>The best practice is simple: use a javadoc style comment above the <code>@Input</code> annotation, with the standard first line summary and subsequent detailed discussion of the meaning of the argument.  These are then picked up by the GATKdocs system and added to the standard walker docs, following the standard structure of GATKDocs <code>@Argument</code> docs.  Below is a best practice documentation example from <code>SelectVariants</code>, which accepts a required variant track and two optional discordance and concordance tracks.</p>
<pre><code class="pre_md">public class SelectVariants extends RodWalker&lt;Integer, Integer&gt; {
   /**
     * Variants from this file are sent through the filtering and modifying routines as directed
     * by the arguments to SelectVariants, and finally are emitted.
     */
    @Input(fullName="variant", shortName = "V", doc="Select variants from this VCF file", required=true)
    public RodBinding&lt;VariantContext&gt; variants;

    /**
     * A site is considered discordant if there exists some sample in eval that has a non-reference genotype
     * and either the site isn't present in this track, the sample isn't present in this track,
     * or the sample is called reference in this track.
     */
    @Input(fullName="discordance", shortName = "disc", doc="Output variants that were not called in this Feature comparison track", required=false)
    private RodBinding&lt;VariantContext&gt; discordanceTrack;

    /**
     * A site is considered concordant if (1) we are not looking for specific samples and there is a variant called
     * in both variants and concordance tracks or (2) every sample present in eval is present in the concordance
     * track and they have the sample genotype call.
     */
    @Input(fullName="concordance", shortName = "conc", doc="Output variants that were also called in this Feature comparison track", required=false)
    private RodBinding&lt;VariantContext&gt; concordanceTrack;
}</code class="pre_md"></pre>
<p>Note how much better the above version is compared to the old pre-<code>Rodbinding</code> syntax (code below).  Below you have a required argument variant that doesn't show up as a formal argument in the GATK, different from the conceptually similar <code>@Arguments</code> for <code>discordanceRodName</code> and <code>concordanceRodName</code>, which have no type restrictions.  There's no place to document the variant argument as well, so the system is effectively blind to this essential argument.</p>
<pre><code class="pre_md">@Requires(value={},referenceMetaData=@RMD(name="variant", type=VariantContext.class))
public class SelectVariants extends RodWalker&lt;Integer, Integer&gt; {
    @Argument(fullName="discordance", shortName =  "disc", doc="Output variants that were not called on a ROD comparison track. Use -disc ROD_NAME", required=false)
    private String discordanceRodName = "";

    @Argument(fullName="concordance", shortName =  "conc", doc="Output variants that were also called on a ROD comparison track. Use -conc ROD_NAME", required=false)
    private String concordanceRodName = "";
}</code class="pre_md"></pre>
<h4>RodBinding examples</h4>
<p>In these examples, we have declared two <code>RodBindings</code> in the Walker</p>
<pre><code class="pre_md">@Input(fullName="mask", doc="Input ROD mask", required=false)
public RodBinding&lt;Feature&gt; mask = RodBinding.makeUnbound(Feature.class);

@Input(fullName="comp", doc="Comparison track", required=false)
public List&lt;RodBinding&lt;VariantContext&gt;&gt; comps = new ArrayList&lt;VariantContext&gt;();</code class="pre_md"></pre>
<ul>
<li>
<p>Get the first value</p>
<p><code>Feature f = tracker.getFirstValue(mask)</code></p>
</li>
<li>
<p>Get all of the values at a location</p>
<p><code>Collection&lt;Feature&gt; fs = tracker.getValues(mask, thisGenomeLoc)</code></p>
</li>
<li>
<p>Get all of the features here, regardless of track </p>
<p><code>Collection&lt;Feature&gt; fs = tracker.getValues(Feature.class)</code></p>
</li>
<li>
<p>Determining if an optional RodBinding was provided
.
if ( mask.isBound() )  // writes out the mask header line, if one was provided
hInfo.add(new VCFFilterHeaderLine(MASK_NAME, &quot;Overlaps a user-input mask&quot;));</p>
<p>if ( ! comps.isEmpty() )
logger.info(&quot;At least one comp was provided&quot;)</p>
</li>
</ul>
<h4>Example usage in Queue scripts</h4>
<p>In <a href="http://gatkforums.broadinstitute.org/discussion/1307/queue-pipeline-scripts-qscripts">QScripts</a> when you need to tag a file use the class <code>TaggedFile</code> which extends from <code>java.io.File</code>.</p>
<table class="table table-striped">
<thead>
<tr>
<th style="text-align: left;">Example</th>
<th style="text-align: left;">in the QScript</th>
<th style="text-align: left;">on the Command Line</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Untagged VCF</td>
<td style="text-align: left;"><code>myWalker.variant = new File("my.vcf")</code></td>
<td style="text-align: left;"><code>-V my.vcf</code></td>
</tr>
<tr>
<td style="text-align: left;">Tagged VCF</td>
<td style="text-align: left;"><code>myWalker.variant = new TaggedFile("my.vcf", "VCF")</code></td>
<td style="text-align: left;"><code>-V:VCF my.vcf</code></td>
</tr>
<tr>
<td style="text-align: left;">Tagged VCF</td>
<td style="text-align: left;"><code>myWalker.variant = new TaggedFile("my.vcf", "VCF,custom=value")</code></td>
<td style="text-align: left;"><code>-V:VCF,custom=value my.vcf</code></td>
</tr>
<tr>
<td style="text-align: left;">Labeling a tumor</td>
<td style="text-align: left;"><code>myWalker.input_file :+= new TaggedFile("mytumor.bam", "tumor")</code></td>
<td style="text-align: left;"><code>-I:tumor mytumor.bam</code></td>
</tr>
</tbody>
</table>
<h4>Notes</h4>
<p>No longer need to (or can) use <code>@Requires</code> and <code>@Allows</code> for ROD data.  This system is now retired.</p>