## Bait bias

http://gatkforums.broadinstitute.org/gatk/discussion/6333/bait-bias

<p>Bait bias (single bait bias or reference bias artifact) is a type of artifact that affects data generated through <a href="http://gatkforums.broadinstitute.org/gatk/discussion/6331">hybrid selection</a> methods. </p>
<p>These artifacts occur during or after the target selection step, and correlate with substitution rates that are biased or higher for sites having one base on the reference/positive strand relative to sites having the complementary base on that strand.  For example, a G&gt;T artifact during the target selection step might result in a higher (G&gt;T)/(C&gt;A) substitution rate at sites with a G on the positive strand (and C on the negative), relative to sites with the flip (C positive)/(G negative). This is known as the <strong>&quot;G-Ref&quot;</strong> artifact.</p>