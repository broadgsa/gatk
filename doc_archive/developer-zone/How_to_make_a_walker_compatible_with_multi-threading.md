## How to make a walker compatible with multi-threading

http://gatkforums.broadinstitute.org/gatk/discussion/2867/how-to-make-a-walker-compatible-with-multi-threading

<p>This document provides an overview of what are the steps required to make a walker multi-threadable using the <code>-nct</code> and the <code>-nt</code> arguments, which make use of the <code>NanoSchedulable</code> and <code>TreeReducible</code> interfaces, respectively.</p>
<hr />
<h3>NanoSchedulable / <code>-nct</code></h3>
<p>Providing <code>-nct</code> support requires that you certify that your walker's <code>map()</code> method is thread-safe -- eg., if any data structures are shared across <code>map()</code> calls, access to these must be properly synchronized. Once your <code>map()</code> method is thread-safe, you can implement the <code>NanoSchedulable</code> interface, an empty interface with no methods that just marks your walker as having a <code>map()</code> method that's safe to parallelize:</p>
<pre><code class="pre_md">/**
 * Root parallelism interface.  Walkers that implement this
 * declare that their map function is thread-safe and so multiple
 * map calls can be run in parallel in the same JVM instance.
 */
public interface NanoSchedulable {
}</code class="pre_md"></pre>
<hr />
<h3>TreeReducible / <code>-nt</code></h3>
<p>Providing <code>-nt</code> support requires that both <code>map()</code> and <code>reduce()</code> be thread-safe, and you also need to implement the <code>TreeReducible</code> interface. Implementing <code>TreeReducible</code> requires you to write a <code>treeReduce()</code> method that tells the engine how to combine the results of multiple <code>reduce()</code> calls:</p>
<pre><code class="pre_md">public interface TreeReducible&lt;ReduceType&gt; {
    /**
     * A composite, 'reduce of reduces' function.
     * @param lhs 'left-most' portion of data in the composite reduce.
     * @param rhs 'right-most' portion of data in the composite reduce.
     * @return The composite reduce type.
     */
    ReduceType treeReduce(ReduceType lhs, ReduceType rhs);
}</code class="pre_md"></pre>
<p>This method differs from <code>reduce()</code> in that while <code>reduce()</code> adds the result of a <em>single</em> <code>map()</code> call onto a running total, <code>treeReduce()</code> takes the aggregated results from multiple map/reduce tasks that have been run in parallel and combines them. So, <code>lhs</code> and <code>rhs</code> might each represent the final result from several hundred map/reduce calls.</p>
<p>Example <code>treeReduce()</code> implementation from the UnifiedGenotyper:</p>
<pre><code class="pre_md">public UGStatistics treeReduce(UGStatistics lhs, UGStatistics rhs) {
    lhs.nBasesCallable += rhs.nBasesCallable;
    lhs.nBasesCalledConfidently += rhs.nBasesCalledConfidently;
    lhs.nBasesVisited += rhs.nBasesVisited;
    lhs.nCallsMade += rhs.nCallsMade;
    return lhs;
}</code class="pre_md"></pre>