## Why are some of the annotation values different with VariantAnnotator compared to UG or HC?

http://gatkforums.broadinstitute.org/gatk/discussion/1550/why-are-some-of-the-annotation-values-different-with-variantannotator-compared-to-ug-or-hc

<p>As featured in <a href="http://gatkforums.broadinstitute.org/discussion/1549/variant-annotator-annotations">this forum question</a>.</p>
<p>Two main things account for these kinds of differences, both linked to default behaviors of the tools:</p>
<ul>
<li>
<p>The tools downsample to different depths of coverage</p>
</li>
<li>The tools apply different read filters</li>
</ul>
<p>In both cases, you can end up looking at different sets or numbers of reads, which causes some of the annotation values to be different. It's usually not a cause for alarm. Remember that many of these annotations should be interpreted <em>relatively</em>, not <em>absolutely</em>.</p>