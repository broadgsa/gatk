## Lane, Library, Sample and Cohort -- what do they mean and why are they important?

http://gatkforums.broadinstitute.org/gatk/discussion/3059/lane-library-sample-and-cohort-what-do-they-mean-and-why-are-they-important

<p>There are four major organizational units for next-generation DNA sequencing processes that used throughout the GATK documentation:</p>
<ul>
<li>
<p><strong>Lane:</strong> The basic machine unit for sequencing. The lane reflects the basic independent run of an NGS machine. For Illumina machines, this is the physical sequencing lane. </p>
</li>
<li>
<p><strong>Library:</strong> A unit of DNA preparation that at some point is physically pooled together. Multiple lanes can be run from aliquots from the same library. The DNA library and its preparation is the natural unit that is being sequenced. For example, if the library has limited complexity, then many sequences are duplicated and will result in a high duplication rate across lanes.</p>
</li>
<li>
<p><strong>Sample:</strong> A single individual, such as human CEPH NA12878. Multiple libraries with different properties can be constructed from the original sample DNA source. Throughout our documentation, we treat samples as independent individuals whose genome sequence we are attempting to determine. Note that from this perspective, tumor / normal samples are different despite coming from the same individual.</p>
</li>
<li><strong>Cohort:</strong> A collection of samples being analyzed together. This organizational unit is the most subjective and depends very specifically on the design goals of the sequencing project. For population discovery projects like the 1000 Genomes, the analysis cohort is the ~100 individual in each population. For exome projects with many deeply sequenced samples (e.g., ESP with 800 EOMI samples) we divide up the complete set of samples into cohorts of ~50 individuals for multi-sample analyses.</li>
</ul>
<p>Note that many GATK commands can be run at the lane level, but will give better results seeing all of the data for a single sample, or even all of the data for all samples. Unfortunately, there's a trade-off in computational cost, since running these commands across all of your data simultaneously requires much more computing power. Please see the documentation for each step to understand what is the best way to group or partition your data for that particular process.</p>