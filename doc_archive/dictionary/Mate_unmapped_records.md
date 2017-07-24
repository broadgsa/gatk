## Mate unmapped records

http://gatkforums.broadinstitute.org/gatk/discussion/6976/mate-unmapped-records

<h3>Mate unmapped records are identifiable using the <code>8</code> SAM flag.</h3>
<p>It is possible for a BAM to have multiple types of mate-unmapped records. These mate unmapped records are distinct from mate missing records, where the mate is altogether absent from the BAM. Of the three types of mate unmapped records listed below, we describe only the first two in this dictionary entry.</p>
<ol>
<li>Singly mapping pair. </li>
<li>A secondary/supplementary record is flagged as mate-unmapped but the mate is in fact mapped.</li>
<li>Both reads in a pair are unmapped.</li>
</ol>
<hr />
<h3>(1) Singly mapping pair</h3>
<p>A mapped read's unmapped mate is marked in their SAM record in an unexpected manner that allow the pair to sort together. If you look at these unmapped reads, the alignment columns 2 and 3 indicate they align, in fact identically to the mapped mate. However, what is distinct is the asterisk <code>*</code> in the CIGAR field (column 6) that indicates the record is unmapped. This allows us to (i) identify the unmapped read as having passed through the aligner, and (ii) keep the pairs together in file manipulations that use either coordinate or queryname sorted BAMs. For example, when a genomic interval of reads are taken to create a new BAM, the pair remain together. For file manipulations dependent on such sorting, we can deduce that these mate unmapped records are immune to becoming missing mates.</p>
<h3>(2) Mate unmapped record whose mate is mapped but in a pair that excludes the record</h3>
<p>The second type of mate unmapped records apply to multimapping read sets processed through MergeBamAlignment such as in <a href="http://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently#latest">Tutorial#6483</a>. Besides reassigning primary and secondary flags within multimapping sets according to a user specified strategy, MergeBamAlignment marks secondary records with the mate unmapped flag. Specifically, after BWA-MEM alignment, records in multimapping sets are all each <em>mate-mapped</em>. After going through MergeBamAlignment, the secondary records become <em>mate-unmapped</em>. The primary alignments remain <em>mate-mapped</em>. This effectively minimizes the association between secondary records from their previous mate. </p>
<hr />
<h3>How do tools treat them differently?</h3>
<p>GATK tools typically ignore secondary/supplementary records from consideration. However, tools will process the mapped read in a singly mapping pair. For example, MarkDuplicates skips secondary records from consideration but marks duplicate singly mapping reads.</p>