## Paired-end / mate-pair

http://gatkforums.broadinstitute.org/gatk/discussion/6327/paired-end-mate-pair

<p>In paired-end sequencing, the library preparation yields a set of fragments, and the machine sequences each fragment from both ends; for example if you have a 300bp contiguous fragment, the machine will sequence e.g. bases 1-75 (forward direction) and bases 225-300 (reverse direction) of the fragment.  </p>
<p>In mate-pair sequencing, the library preparation yields two fragments that are distal to each other in the genome and in the opposite in orientation to that of a mate-paired fragment.</p>
<p>The three read orientation categories are forward reverse (FR), reverse forward (RF), and reverse-reverse/forward-forward (TANDEM). In general, paired-end reads tend to be in a FR orientation, have relatively small inserts (~300 - 500 bp), and are particularly useful for the sequencing of fragments that contain short repeat regions.  Mate-pair fragments are generally in a RF conformation, contain larger inserts (~3 kb), and enable sequence coverage of genomic regions containing large structural rearrangements. Tandem reads can result from inversions and rearrangements during library preparation. </p>
<p>Here is a more illustrative example:</p>
<p><strong>FR:</strong> 5' --F--&gt;       &lt;--R-- 5' (in slang called &quot;innie&quot; because they point inward)</p>
<p><strong>RF:</strong> &lt;--R-- 5'       5' --F--&gt; (in slang called &quot;outie&quot; because they point outward)</p>
<p><strong>TANDEM:</strong> 5' --F--&gt;   5' --F--&gt;  or  &lt;--R-- 5'   &lt;--R-- 5'</p>
<p>The figure below illustrates this graphically along with the SAM flags that correspond to the FR and RF configurations.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/e3/c9e87118d6e8c4b2a4e014d97a1b22.png" />
<p>For detailed explanations of library construction strategies (for Illumina sequencers) and how read orientations are determined, please see:</p>
<ul>
<li><a href="http://www.illumina.com/technology/next-generation-sequencing/paired-end-sequencing_assay.html">Illumina paired-end sequencing documentation (webpage)</a></li>
<li><a href="http://www.illumina.com/documents/products/technotes/technote_nextera_matepair_data_processing.pdf">Illumina Nextera mate-pair processing documentation (pdf)</a></li>
</ul>