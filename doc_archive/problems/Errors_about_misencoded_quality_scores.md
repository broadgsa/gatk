## Errors about misencoded quality scores

http://gatkforums.broadinstitute.org/gatk/discussion/6470/errors-about-misencoded-quality-scores

<h3>The problem</h3>
<p>You get an error like this: </p>
<pre><code class="pre_md">SAM/BAM/CRAM file &lt;filename&gt; appears to be using the wrong encoding for quality scores</code class="pre_md"></pre>
<h3>Why this happens</h3>
<p>The standard format for quality score encodings is that Q0 == ASCII 33 according to the SAM specification. However, in some datasets (including older Illumina data), encoding starts at ASCII 64. This is a problem because the GATK assumes that it can use the quality scores as they are. If they are in fact encoded using a different scale, our tools will make an incorrect estimation of the quality of your data, and your analysis results will be off. </p>
<p>To prevent this from happening, the GATK engine performs a sanity check of the quality score encodings that will abort the program run if they are not standard (since version 2.3), and output the error message shown above. </p>
<h3>Solution</h3>
<p>If this happens to you, you'll need to run again with the flag <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_engine_CommandLineGATK.php#--fix_misencoded_quality_scores"> <code>--fix_misencoded_quality_scores</code> / <code>-fixMisencodedQuals</code></a>. What will happen is that the engine will simply subtract 31 from every quality score as it is read in, and proceed with the corrected values. Output files will include the correct scores where applicable.</p>
<h3>Related problems</h3>
<p>In some cases the data contains a mix of encodings (which is likely to arise if you're passing in a lot of different files from different sources together), and the GATK can't automatically compensate for that. There is an argument you can use to override this check: <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_engine_CommandLineGATK.php#"><code>-allowPotentiallyMisencodedQuals</code> / <code>--allow_potentially_misencoded_quality_scores</code></a>; but you use it at your own risk. We strongly encourage you to check the encodings of your files rather than use this option. </p>