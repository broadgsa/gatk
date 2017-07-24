## Writing GATKdocs for your walkers

http://gatkforums.broadinstitute.org/gatk/discussion/1324/writing-gatkdocs-for-your-walkers

<p>The GATKDocs are what we call <a href="http://www.broadinstitute.org/gatk/gatkdocs/">&quot;Technical Documentation&quot;</a> in the Guide section of this website. The HTML pages are generated automatically at build time from specific blocks of documentation in the source code. </p>
<p>The best place to look for example documentation for a GATK walker is GATKDocsExample walker in <code>org.broadinstitute.sting.gatk.examples</code>.  This is available <a href="https://github.com/broadgsa/gatk/blob/master/public/java/src/org/broadinstitute/sting/gatk/examples/GATKDocsExample.java">here</a>.  </p>
<p>Below is the reproduction of that file from August 11, 2011:</p>
<pre><code class="pre_md">/**
 * [Short one sentence description of this walker]
 *
 * &lt;p&gt;
 * [Functionality of this walker]
 * &lt;/p&gt;
 *
 * &lt;h2&gt;Input&lt;/h2&gt;
 * &lt;p&gt;
 * [Input description]
 * &lt;/p&gt;
 *
 * &lt;h2&gt;Output&lt;/h2&gt;
 * &lt;p&gt;
 * [Output description]
 * &lt;/p&gt;
 *
 * &lt;h2&gt;Examples&lt;/h2&gt;
 * PRE-TAG
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T $WalkerName
 * PRE-TAG
 *
 * @category Walker Category
 * @author Your Name
 * @since Date created
 */
public class GATKDocsExample extends RodWalker&lt;Integer, Integer&gt; {
    /**
     * Put detailed documentation about the argument here.  No need to duplicate the summary information
     * in doc annotation field, as that will be added before this text in the documentation page.
     *
     * Notes:
     * &lt;ul&gt;
     *     &lt;li&gt;This field can contain HTML as a normal javadoc&lt;/li&gt;
     *     &lt;li&gt;Don't include information about the default value, as gatkdocs adds this automatically&lt;/li&gt;
     *     &lt;li&gt;Try your best to describe in detail the behavior of the argument, as ultimately confusing
     *          docs here will just result in user posts on the forum&lt;/li&gt;
     * &lt;/ul&gt;
     */
    @Argument(fullName="full", shortName="short", doc="Brief summary of argument [~ 80 characters of text]", required=false)
    private boolean myWalkerArgument = false;

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) { return 0; }
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer value, Integer sum) { return value + sum; }
    public void onTraversalDone(Integer result) { }
}</code class="pre_md"></pre>