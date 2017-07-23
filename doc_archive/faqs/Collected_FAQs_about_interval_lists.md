## Collected FAQs about interval lists

http://gatkforums.broadinstitute.org/gatk/discussion/1319/collected-faqs-about-interval-lists

<h2>1. Can GATK tools be restricted to specific intervals instead of processing the entire reference?</h2>
<p>Absolutely. Just use the <code>-L</code> argument to provide the list of intervals you wish to run on. Or you can use <code>-XL</code> to <em>exclude</em> intervals, e.g. to blacklist genome regions that are problematic. </p>
<hr />
<h2>2. What file formats does GATK support for interval lists?</h2>
<p>GATK supports several types of interval list formats: Picard-style <code>.interval_list</code>, GATK-style <code>.list</code>, BED files with extension <code>.bed</code>, and VCF files.  </p>
<h3>A. Picard-style <code>.interval_list</code></h3>
<p>Picard-style interval files have a SAM-like header that includes a sequence dictionary. The intervals are given in the form <code>&lt;chr&gt; &lt;start&gt; &lt;stop&gt; + &lt;target_name&gt;</code>, with fields separated by tabs, and the coordinates are 1-based (first position in the genome is position 1, not position 0). </p>
<pre><code class="pre_md">@HD     VN:1.0  SO:coordinate
@SQ     SN:1    LN:249250621    AS:GRCh37       UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta   M5:1b22b98cdeb4a9304cb5d48026a85128     SP:Homo Sapiens
@SQ     SN:2    LN:243199373    AS:GRCh37       UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta   M5:a0d9851da00400dec1098a9255ac712e     SP:Homo Sapiens
1       30366   30503   +       target_1
1       69089   70010   +       target_2
1       367657  368599  +       target_3
1       621094  622036  +       target_4
1       861320  861395  +       target_5
1       865533  865718  +       target_6</code class="pre_md"></pre>
<p>This is the preferred format because the explicit sequence dictionary safeguards against accidental misuse (e.g. apply hg18 intervals to an hg19 BAM file). Note that this file is 1-based, not 0-based (the first position in the genome is position 1).</p>
<h3>B. GATK-style <code>.list</code> or <code>.intervals</code></h3>
<p>This is a simpler format, where intervals are in the form <code>&lt;chr&gt;:&lt;start&gt;-&lt;stop&gt;</code>, and no sequence dictionary is necessary. This file format also uses 1-based coordinates. Note that only the <code>&lt;chr&gt;</code> part is strictly required; if you just want to specify chromosomes/ contigs as opposed to specific coordinate ranges, you don't need to specify the rest. Both <code>&lt;chr&gt;:&lt;start&gt;-&lt;stop&gt;</code> and <code>&lt;chr&gt;</code> can be present in the same file. You can also specify intervals in this format directly at the command line instead of writing them in a file.</p>
<h3>C. BED files with extension <code>.bed</code></h3>
<p>We also accept the widely-used BED format, where intervals are in the form <code>&lt;chr&gt; &lt;start&gt; &lt;stop&gt;</code>, with fields separated by tabs. However, you should be aware that this file format is 0-based for the start coordinates, so coordinates taken from 1-based formats (e.g. if you're cooking up a custom interval list derived from a file in a 1-based format) should be offset by 1. The GATK engine recognizes the <code>.bed</code> extension and interprets the coordinate system accordingly.</p>
<h3>D. VCF files</h3>
<p>Yeah, I bet you didn't expect that was a thing! It's very convenient. Say you want to redo a variant calling run on a set of variant calls that you were given by a colleague, but with the latest version of HaplotypeCaller. You just provide the VCF, slap on some padding on the fly using e.g. <code>-ip 100</code> in the HC command, and boom, done. Each record in the VCF will be interpreted as a single-base interval, and by adding padding you ensure that the caller sees enough context to reevaluate the call appropriately.</p>
<hr />
<h2>3. Is there a required order of intervals?</h2>
<p>Yes, thanks for asking. The intervals MUST be sorted by coordinate (in increasing order) within contigs; and the contigs must be sorted in the same order as in the sequence dictionary. This is for efficiency reasons. </p>
<hr />
<h2>4. Can I provide multiple sets of intervals?</h2>
<p>Sure, no problem -- just pass them in using separate <code>-L</code> arguments. You can use all the different formats within the same command line. By default, the GATK engine will take the UNION of all the intervals in all the sets. This behavior can be modified by setting an <a href="https://software.broadinstitute.org/gatk/documentation/tooldocs/org_broadinstitute_gatk_engine_CommandLineGATK.php#--interval_set_rule"><code>interval_set</code></a> rule.</p>
<hr />
<h2>5. How will GATK handle intervals that abut or overlap?</h2>
<p>Very gracefully. By default the GATK engine will merge any intervals that abut (i.e. they are contiguous, they touch without overlapping) or overlap into a single interval. This behavior can be modified by setting an <a href="https://software.broadinstitute.org/gatk/documentation/tooldocs/org_broadinstitute_gatk_engine_CommandLineGATK.php#--interval_merging"><code>interval_merging</code></a> rule.</p>
<hr />
<h2>6. What's the best way to pad intervals?</h2>
<p>You can use the <a href="https://software.broadinstitute.org/gatk/documentation/tooldocs/org_broadinstitute_gatk_engine_CommandLineGATK.php#--interval_padding"><code>-ip</code></a> engine argument to add padding on the fly. No need to produce separate padded targets files. Sweet, right? </p>
<p>Note that if intervals that previously didn't abut or overlap before you added padding now do so, by default the GATK engine will merge them as described above. This behavior can be modified by setting an <a href="https://software.broadinstitute.org/gatk/documentation/tooldocs/org_broadinstitute_gatk_engine_CommandLineGATK.php#--interval_merging"><code>interval_merging</code></a> rule.</p>