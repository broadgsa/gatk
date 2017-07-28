## Local Realignment around Indels

http://gatkforums.broadinstitute.org/gatk/discussion/38/local-realignment-around-indels

<h4>For a discussion of the implications of removing indel realignment from workflows, see <a href="http://gatkforums.broadinstitute.org/gatk/discussion/7847">Blog#7847</a> from June 2016.</h4>
<hr />
<h2>Realigner Target Creator </h2>
<p>For a complete, detailed argument reference, refer to the GATK document page <a rel="nofollow" class="external text" href="http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_indels_RealignerTargetCreator.html">here</a>.
</p><p><br />
</p>
<h2>Indel Realigner</h2>
<p>For a complete, detailed argument reference, refer to the GATK document page <a rel="nofollow" class="external text" href="http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_indels_IndelRealigner.html">here</a>.
</p><p><br />
</p>
<hr />
<h1>Running the Indel Realigner only at known sites</h1>
<p>While we advocate for using the Indel Realigner over an aggregated bam using the full Smith-Waterman alignment algorithm, it will work for just a single lane of sequencing data when run in -knownsOnly mode.  Novel sites obviously won't be cleaned up, but the majority of a single individual's short indels will already have been seen in dbSNP and/or 1000 Genomes.  One would employ the known-only/lane-level realignment strategy in a large-scale project (e.g. 1000 Genomes) where computation time is severely constrained and limited.  We modify the example arguments from above to reflect the command-lines necessary for known-only/lane-level cleaning.
</p><p>The RealignerTargetCreator step would need to be done just once for a single set of indels; so as long as the set of known indels doesn't change, the output.intervals file from below would never need to be recalculated.
</p>
<pre>
 java -Xmx1g -jar /path/to/GenomeAnalysisTK.jar \
  -T RealignerTargetCreator \
  -R /path/to/reference.fasta \
  -o /path/to/output.intervals \
  -known /path/to/indel_calls.vcf
</pre>
<p>The IndelRealigner step needs to be run on every bam file.
</p>
<pre>
java -Xmx4g -Djava.io.tmpdir=/path/to/tmpdir \
  -jar /path/to/GenomeAnalysisTK.jar \
  -I &lt;lane-level.bam&gt; \
  -R &lt;ref.fasta&gt; \
  -T IndelRealigner \
  -targetIntervals &lt;intervalListFromStep1Above.intervals&gt; \
  -o &lt;realignedBam.bam&gt; \
  -known /path/to/indel_calls.vcf
  --consensusDeterminationModel KNOWNS_ONLY \
  -LOD 0.4
</pre>