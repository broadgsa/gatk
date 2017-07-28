## PF reads / Illumina chastity filter

http://gatkforums.broadinstitute.org/gatk/discussion/6329/pf-reads-illumina-chastity-filter

<p>Illumina sequencers perform an internal quality filtering procedure called <strong>chastity filter</strong>, and reads that pass this filter are called <strong>PF</strong> for <strong>pass-filter</strong>. According to Illumina, <strong>chastity</strong> is defined as the ratio of the brightest base intensity divided by the sum of the brightest and second brightest base intensities. Clusters of reads pass the filter if no more than 1 base call has a chastity value below 0.6 in the first 25 cycles. This filtration process removes the least reliable clusters from the image analysis results.</p>
<p>For additional information on chastity filters, please see:  </p>
<ul>
<li>Illumina, Inc. (2015).  Calculating Percent Passing Filter for Patterned and Non-Patterned Flow Cells: A comparison of methods for calculating percent passing filter on Illumina flow cells</li>
<li>Ilumina Inc. (2014) HiSeq X System user guide</li>
</ul>
<p>Both articles can be found at <a href="http://www.Illumina.com">http://www.Illumina.com</a></p>