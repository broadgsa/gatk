package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import net.sf.samtools.SAMRecord;

import java.util.List;

public class Pilot3CoverageWalker extends LocusWalker<Integer, Integer> {

    @Argument(fullName = "extended", shortName="ext", doc="extended output", required=false)
    public boolean extendedOutput = false;

    @Argument(fullName="min_mapq", shortName="mmq", required=false, doc="Minimum mapping quality of reads to consider") public Integer MIN_MAPQ = 1;


    public void initialize() {
        out.println("track type=wiggle_0 name=Pilot3Coverage viewLimits=0:1");
    }

    public String walkerType() { return "ByLocus"; }

    // Do we actually want to operate on the context?
    public boolean filter(RefMetaDataTracker tracker, char ref, AlignmentContext context) {
        return true;    // We are keeping all the reads
    }



    // Map over the org.broadinstitute.sting.gatk.contexts.AlignmentContext
    int MAPPING_QUALITY_THRESHOLD = 1;
    int totalSites;
    int tumorCovered;
    int normalCovered;
    int somaticCovered;
    long start = 0;

    int lastContigIndex = -1;
    long lastPosition = -1;

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (start ==0) { start = System.currentTimeMillis(); }

        List<SAMRecord> reads = context.getReads();
        int totalDepth = 0;
        int positiveStrandDepth = 0;
        int negativeStrandDepth = 0;
        for ( int i = 0; i < reads.size(); i++ )
        {
            SAMRecord read = reads.get(i);

            // TODO: could this be done better elsewhere?
            // only process primary, non duplicate alignments
            // that come from fully mapped pairs with a mappign quality threshold >= x
            if (read.getNotPrimaryAlignmentFlag() ||
                read.getDuplicateReadFlag() ||
                read.getReadUnmappedFlag() ||
                read.getMappingQuality() < MIN_MAPQ
//        ||
//                read.getMateUnmappedFlag() ||
                    ) {
                continue;
            }



            totalDepth++;
            if (read.getReadNegativeStrandFlag()) {
                negativeStrandDepth++;
            } else {
                positiveStrandDepth++;
            }

        }


        // if the contig index has changed OR if it's the same contig index but we jumped positions
        // output a wiggle header
        StringBuilder sb = new StringBuilder();
        if (lastContigIndex != context.getLocation().getContigIndex() ||
            lastPosition + 1 != context.getPosition()) {
                lastContigIndex = context.getLocation().getContigIndex();
                sb.append("fixedStep").append(" ")
                  .append("chrom=").append(context.getContig()).append(" ")
                  .append("start=").append(context.getPosition()).append(" ")
                  .append("step=1")
                  .append("\n");
        }
        lastPosition = context.getPosition();

        if (extendedOutput) {
            sb.append(context.getPosition()).append(" ");
            sb.append(totalDepth).append(" ");
            sb.append(positiveStrandDepth).append(" ");
            sb.append(negativeStrandDepth).append(" ");
        }

        boolean siteCovered = (totalDepth >= 10 && positiveStrandDepth >= 2 && negativeStrandDepth >= 2);
        sb.append((siteCovered)?"1":"0");

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
