package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.HapMapAlleleFrequenciesROD;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import net.sf.samtools.SAMRecord;

import java.util.List;
import java.util.Formatter;
import static java.lang.Math.log10;

import edu.mit.broad.picard.genotype.DiploidGenotype;
import edu.mit.broad.picard.PicardException;

public class SomaticCoverageWalker extends LocusWalker<Integer, Integer> {

    @Argument(fullName = "tumor_sample_name", shortName = "1", required = true)
    String tumorSampleName = "TCGA-06-0188-01A-01W";  // FIXME: cmdline parsing not working!

    @Argument(fullName = "normal_sample_name", shortName = "2", required = true)
    String normalSampleName = "TCGA-06-0188-10B-01W"; // FIXME: cmdline parsing not working!

    public void initialize() {
    }

    public String walkerType() { return "ByLocus"; }

    // Do we actually want to operate on the context?
    public boolean filter(RefMetaDataTracker tracker, char ref, LocusContext context) {
        return true;    // We are keeping all the reads
    }



    // Map over the org.broadinstitute.sting.gatk.LocusContext
    int MAPPING_QUALITY_THRESHOLD = 1;
    int totalSites;
    int tumorCovered;
    int normalCovered;
    int somaticCovered;
    long start = 0;

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        if (start ==0) { start = System.currentTimeMillis(); }

        List<SAMRecord> reads = context.getReads();
        int tumorDepth = 0;
        int normalDepth = 0;
        for ( int i = 0; i < reads.size(); i++ )
        {
            SAMRecord read = reads.get(i);

            // TODO: could this be done better elsewhere?
            // only process primary, non duplicate alignments
            // that come from fully mapped pairs with a mappign quality threshold >= x
            if (read.getNotPrimaryAlignmentFlag() ||
                read.getDuplicateReadFlag() ||
                read.getReadUnmappedFlag() ||
                read.getMappingQuality() < MAPPING_QUALITY_THRESHOLD
//        ||
//                read.getMateUnmappedFlag() ||
                    ) {
                continue;
            }

            String rg = (String) read.getAttribute("RG");
            String sample = read.getHeader().getReadGroup(rg).getSample();



            if (normalSampleName.equals(sample)) {
                normalDepth++;
            } else if (tumorSampleName.equals(sample)) {
                tumorDepth++;
            } else {
                throw new RuntimeException("Unknown Sample Name: " + sample);
            }

        }


        boolean isTumorCovered = tumorDepth >= 14;
        boolean isNormalCovered = normalDepth >= 8;

        if (isTumorCovered) { tumorCovered++; }
        if (isNormalCovered) { normalCovered++; }
        if (isTumorCovered && isNormalCovered) {somaticCovered++; }
        totalSites++;

//        if (totalSites % 100000 == 0) {
//            long now = System.currentTimeMillis();
//            out.println(String.format("%s:%d %d %d %d %d %dms", context.getContig(), context.getPosition(), totalSites, tumorCovered, normalCovered, somaticCovered, (now-start)));
//        }

        StringBuilder sb = new StringBuilder();
        sb.append(context.getContig()).append(" ");
        sb.append(context.getPosition());
        out.println(sb.toString());
        return 1;
    }

    // Given result of map function
    public Integer reduceInit() {
        return 0;
    }
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    @Override
    public void onTraversalDone(Integer result) {
//        out.println(String.format("FINAL - %d %d %d %d", totalSites, tumorCovered, normalCovered, somaticCovered));
    }



}