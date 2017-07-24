## HC step 1: Defining ActiveRegions by measuring data entropy

http://gatkforums.broadinstitute.org/gatk/discussion/4147/hc-step-1-defining-activeregions-by-measuring-data-entropy

<p>This document describes the procedure used by HaplotypeCaller to define ActiveRegions on which to operate as a prelude to variant calling. For more context information on how this fits into the overall HaplotypeCaller method, please see the more general <a href="http://www.broadinstitute.org/gatk/guide/article?id=4148">HaplotypeCaller documentation</a>.</p>
<h3>Summary</h3>
<p>To define active regions, the HaplotypeCaller operates in three phases. First, it computes an <strong>activity score</strong> for each individual genome position, yielding the <strong>raw activity profile</strong>, which is a wave function of activity per position. Then, it applies a smoothing algorithm to the raw profile, which is essentially a sort of averaging process, to yield the actual <strong>activity profile</strong>. Finally, it identifies local maxima where the activity profile curve rises above the preset activity threshold, and defines appropriate intervals to encompass the active profile within the preset size constraints. </p>
<hr />
<h3>1. Calculating the raw activity profile</h3>
<p>Active regions are determined by calculating a profile function that characterizes “interesting” regions likely to contain variants. The raw profile is first calculated locus by locus. </p>
<p>In the normal case (no special mode is enabled) the per-position score is the probability that the position contains a variant as calculated using the reference-confidence model applied to the original alignment.</p>
<p>If using the mode for genotyping given alleles (GGA) or the advanced-level flag <code>-useAlleleTrigger</code>, and the site is overlapped by an allele in the VCF file provided through the <code>-alleles</code> argument, the score is set to 1. If the position is not covered by a provided allele, the score is set to 0. </p>
<p>This operation gives us a single raw value for each position on the genome (or within the analysis intervals requested using the <code>-L</code> argument).</p>
<hr />
<h3>2. Smoothing the activity profile</h3>
<p>The final profile is calculated by smoothing this initial raw profile following three steps. The first two steps consist in spreading individual position raw profile values to contiguous bases. As a result each position will have more than one raw profile value that are added up in the third and last step to obtain a final unique and smoothed value per position.</p>
<ol>
<li>
<p>Unless one of the special modes is enabled (GGA or allele triggering), the position profile value will be copied over to adjacent regions if enough high quality soft-clipped bases immediately precede or follow that position in the original alignment. At time of writing, high-quality soft-clipped bases are those with quality score of Q29 or more. We consider that there are enough of such a soft-clips when the average number of high quality bases per soft-clip is 7 or more. In this case the site profile value is copied to all bases within a radius of that position as large as the average soft-clip length without exceeding a maximum of 50bp.</p>
</li>
<li>
<p>Each profile value is then divided and spread out using a Gaussian kernel covering up to 50bp radius centered at its current position with a standard deviation, or sigma, set using the <code>-bandPassSigma</code> argument (current default is 17 bp). The larger the sigma, the broader the spread will be.</p>
</li>
<li>For each position, the final smoothed value is calculated as the sum of all its profile values after steps 1 and 2.</li>
</ol>
<hr />
<h3>3. Setting the ActiveRegion thresholds and intervals</h3>
<p>The resulting profile line is cut in regions where it crosses the non-active to active threshold (currently set to 0.002). Then we make some adjustments to these boundaries so that those regions that are to be considered active, with a profile running over that threshold, fall within the minimum (fixed to 50bp) and maximum region size (customizable using <code>-activeRegionMaxSize</code>).</p>
<ul>
<li>
<p>If the region size falls within the limits we leave it untouched (it's good to go).</p>
</li>
<li>
<p>If the region size is shorter than the minimum, it is greedily extended forward ignoring that cut point and we come back to step 1. Only if this is not possible because we hit a hard-limit (end of the chromosome or requested analysis interval) we will accept the small region as it is.</p>
</li>
<li>If it is too long, we find the lowest local minimum between the maximum and minimum region size. A local minimum is a profile value preceded by a large one right up-stream (-1bp) and an equal or larger value down-stream (+1bp). In case of a tie, the one further downstream takes precedence. If there is no local minimum we simply force the cut so that the region has the maximum active region size. </li>
</ul>
<p>Of the resulting regions, those with a profile that runs over this threshold are considered active regions and progress to variant discovery and or calling whereas regions whose profile runs under the threshold are considered inactive regions and are discarded except if we are running HC in reference confidence mode.</p>
<p>There is a final post-processing step to clean up and trim the ActiveRegion:</p>
<ul>
<li>
<p>Remove bases at each end of the read (hard-clipping) until there a base with a call quality equal or greater than minimum base quality score (customizable parameter <code>-mbq</code>, 10 by default).  </p>
</li>
<li>
<p>Include or exclude remaining soft-clipped ends. Soft clipped ends will be used for assembly and calling unless the user has requested their exclusion (using <code>-dontUseSoftClippedBases</code>), if the read and its mate map to the same chromosome, and if they are in the correct standard orientation (i.e. LR and RL).</p>
</li>
<li>
<p>Clip off adaptor sequences of the read if present.</p>
</li>
<li>
<p>Discard all reads that no longer overlap with the ActiveRegion after the trimming operations described above.</p>
</li>
<li>Downsample remaining reads to a maximum of 1000 reads per sample, but respecting a minimum of 5 reads starting per position. This is performed after any downsampling by the traversal itself (<code>-dt</code>, <code>-dfrac</code>, <code>-dcov</code> etc.) and cannot be overriden from the command line.  </li>
</ul>