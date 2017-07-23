## How can I use parallelism to make GATK tools run faster?

http://gatkforums.broadinstitute.org/gatk/discussion/1975/how-can-i-use-parallelism-to-make-gatk-tools-run-faster

<p><em>This document provides technical details and recommendations on how the parallelism options offered by the GATK can be used to yield optimal performance results.</em></p>
<h3>Overview</h3>
<p>As explained in the <a href="http://www.broadinstitute.org/gatk/guide/article?id=1988">primer on parallelism for the GATK</a>, there are two main kinds of parallelism that can be applied to the GATK: multi-threading and scatter-gather (using <a href="https://software.broadinstitute.org/gatk/documentation/pipelines">Queue or Crom/WDL</a>).</p>
<h3>Multi-threading options</h3>
<p>There are two options for multi-threading with the GATK, controlled by the arguments <code>-nt</code> and <code>-nct</code>, respectively, which can be combined:</p>
<ul>
<li><code>-nt / --num_threads</code>
controls the number of <strong>data threads</strong> sent to the processor </li>
<li><code>-nct / --num_cpu_threads_per_data_thread</code>
controls the number of <strong>CPU threads</strong> allocated to each data thread</li>
</ul>
<p>For more information on how these multi-threading options work, please read the <a href="http://www.broadinstitute.org/gatk/guide/article?id=1988">primer on parallelism for the GATK</a>.</p>
<h4>Memory considerations for multi-threading</h4>
<p>Each data thread needs to be given the full amount of memory you’d normally give a single run. So if you’re running a tool that normally requires 2 Gb of memory to run, if you use <code>-nt 4</code>, the multithreaded run will  use 8 Gb of memory. In contrast, CPU threads will share the memory allocated to their “mother” data thread, so you don’t need to worry about allocating memory based on the number of CPU threads you use. </p>
<h4>Additional consideration when using <code>-nct</code> with versions 2.2 and 2.3</h4>
<p>Because of the way the <code>-nct</code> option was originally implemented, in versions 2.2 and 2.3, there is one CPU thread that is reserved by the system to “manage” the rest. So if you use <code>-nct</code>, you’ll only really start seeing a speedup with <code>-nct 3</code> (which yields two effective &quot;working&quot; threads) and above. This limitation has been resolved in the implementation that will be available in versions 2.4 and up.</p>
<h3>Scatter-gather</h3>
<p>For more details on scatter-gather, see the <a href="http://gatkforums.broadinstitute.org/discussion/1988/a-primer-on-parallelism-with-the-gatk">primer on parallelism for the GATK</a> and the documentation on <a href="https://software.broadinstitute.org/gatk/documentation/pipelines">pipelining options</a>.</p>
<h3>Applicability of parallelism to the major GATK tools</h3>
<p>Please note that not all tools support all parallelization modes. The parallelization modes that are available for each tool depend partly on the type of traversal that the tool uses to walk through the data, and partly on the nature of the analyses it performs.</p>
<table class="table table-striped">
<thead>
<tr>
<th style="text-align: left;">Tool</th>
<th style="text-align: left;">Full name</th>
<th style="text-align: left;">Type of traversal</th>
<th style="text-align: center;">NT</th>
<th style="text-align: center;">NCT</th>
<th style="text-align: center;">SG</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">RTC</td>
<td style="text-align: left;">RealignerTargetCreator</td>
<td style="text-align: left;">RodWalker</td>
<td style="text-align: center;">+</td>
<td style="text-align: center;">-</td>
<td style="text-align: center;">-</td>
</tr>
<tr>
<td style="text-align: left;">IR</td>
<td style="text-align: left;">IndelRealigner</td>
<td style="text-align: left;">ReadWalker</td>
<td style="text-align: center;">-</td>
<td style="text-align: center;">-</td>
<td style="text-align: center;">+</td>
</tr>
<tr>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">BaseRecalibrator</td>
<td style="text-align: left;">LocusWalker</td>
<td style="text-align: center;">-</td>
<td style="text-align: center;">+</td>
<td style="text-align: center;">+</td>
</tr>
<tr>
<td style="text-align: left;">PR</td>
<td style="text-align: left;">PrintReads</td>
<td style="text-align: left;">ReadWalker</td>
<td style="text-align: center;">-</td>
<td style="text-align: center;">+</td>
<td style="text-align: center;">-</td>
</tr>
<tr>
<td style="text-align: left;">RR</td>
<td style="text-align: left;">ReduceReads</td>
<td style="text-align: left;">ReadWalker</td>
<td style="text-align: center;">-</td>
<td style="text-align: center;">-</td>
<td style="text-align: center;">+</td>
</tr>
<tr>
<td style="text-align: left;">HC</td>
<td style="text-align: left;">HaplotypeCaller</td>
<td style="text-align: left;">ActiveRegionWalker</td>
<td style="text-align: center;">-</td>
<td style="text-align: center;">(+)</td>
<td style="text-align: center;">+</td>
</tr>
<tr>
<td style="text-align: left;">UG</td>
<td style="text-align: left;">UnifiedGenotyper</td>
<td style="text-align: left;">LocusWalker</td>
<td style="text-align: center;">+</td>
<td style="text-align: center;">+</td>
<td style="text-align: center;">+</td>
</tr>
</tbody>
</table>
<p>Note that while HaplotypeCaller supports <code>-nct</code> in principle, many have reported that it is not very stable (random crashes may occur -- but if there is no crash, results will be correct). We prefer not to use this option with HC; use it at your own risk. </p>
<h3>Recommended configurations</h3>
<p>The table below summarizes configurations that we typically use for our own projects (one per tool, except we give three alternate possibilities for the UnifiedGenotyper). The different values allocated for each tool reflect not only the technical capabilities of these tools (which options are supported), but also our empirical observations of what provides the best tradeoffs between performance gains and commitment of resources. Please note however that this is meant only as a guide, and that we cannot give you any guarantee that these configurations are the best for your own setup. You will probably have to experiment with the settings to find the configuration that is right for you. </p>
<table class="table table-striped">
<thead>
<tr>
<th style="text-align: left;">Tool</th>
<th style="text-align: center;">RTC</th>
<th style="text-align: center;">IR</th>
<th style="text-align: center;">BR</th>
<th style="text-align: center;">PR</th>
<th style="text-align: center;">RR</th>
<th style="text-align: center;">HC</th>
<th style="text-align: center;">UG</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Available modes</td>
<td style="text-align: center;">NT</td>
<td style="text-align: center;">SG</td>
<td style="text-align: center;">NCT,SG</td>
<td style="text-align: center;">NCT</td>
<td style="text-align: center;">SG</td>
<td style="text-align: center;">NCT,SG</td>
<td style="text-align: center;">NT,NCT,SG</td>
</tr>
<tr>
<td style="text-align: left;">Cluster nodes</td>
<td style="text-align: center;">1</td>
<td style="text-align: center;">4</td>
<td style="text-align: center;">4</td>
<td style="text-align: center;">1</td>
<td style="text-align: center;">4</td>
<td style="text-align: center;">4</td>
<td style="text-align: center;">4 / 4 / 4</td>
</tr>
<tr>
<td style="text-align: left;">CPU threads (<code>-nct</code>)</td>
<td style="text-align: center;">1</td>
<td style="text-align: center;">1</td>
<td style="text-align: center;">8</td>
<td style="text-align: center;">4-8</td>
<td style="text-align: center;">1</td>
<td style="text-align: center;">4</td>
<td style="text-align: center;">3 / 6 / 24</td>
</tr>
<tr>
<td style="text-align: left;">Data threads (<code>-nt</code>)</td>
<td style="text-align: center;">24</td>
<td style="text-align: center;">1</td>
<td style="text-align: center;">1</td>
<td style="text-align: center;">1</td>
<td style="text-align: center;">1</td>
<td style="text-align: center;">1</td>
<td style="text-align: center;">8 / 4 / 1</td>
</tr>
<tr>
<td style="text-align: left;">Memory (Gb)</td>
<td style="text-align: center;">48</td>
<td style="text-align: center;">4</td>
<td style="text-align: center;">4</td>
<td style="text-align: center;">4</td>
<td style="text-align: center;">4</td>
<td style="text-align: center;">16</td>
<td style="text-align: center;">32 / 16 / 4</td>
</tr>
</tbody>
</table>
<p>Where NT is data multithreading, NCT is CPU multithreading and SG is scatter-gather using Queue or other data parallelization framework. For more details on scatter-gather, see the <a href="http://www.broadinstitute.org/gatk/guide/article?id=1988">primer on parallelism for the GATK</a> and the documentation on <a href="https://software.broadinstitute.org/gatk/documentation/pipelines">pipelining options</a>.</p>