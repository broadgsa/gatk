package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.ReadBackedPileup;

import java.util.*;

import net.sf.samtools.SAMRecord;

public class PooledCalculationModel extends JointEstimateGenotypeCalculationModel {

    private static final String POOL_SAMPLE_NAME = "POOL";

    protected PooledCalculationModel() {}

    protected int getNSamples(HashMap<String, AlignmentContextBySample> contexts) {
        return POOL_SIZE;
    }
    
    protected HashMap<String, AlignmentContextBySample> createContexts(AlignmentContext context) {
        // for testing purposes, we may want to throw multi-samples at pooled mode,
        // so we can't use the standard splitContextBySample() method here

        AlignmentContextBySample pooledContext = new AlignmentContextBySample(context.getLocation());

        int deletionsInPileup = 0;
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();

        for (int i = 0; i < reads.size(); i++) {

            // check for deletions
            int offset = offsets.get(i);
            if ( offset == -1 ) {
                // are there too many deletions in the pileup?
                if ( ++deletionsInPileup > maxDeletionsInPileup && maxDeletionsInPileup >= 0 )
                    return null;
            }

            // add the read to this sample's context
            // note that bad bases are added to the context (for DoC calculations later)
            pooledContext.add(reads.get(i), offset);
        }

        HashMap<String, AlignmentContextBySample> contexts = new HashMap<String, AlignmentContextBySample>();
        contexts.put(POOL_SAMPLE_NAME, pooledContext);
        return contexts;
    }
    
    protected double computeLog10PofDgivenAFi(char ref, char alt, double f, HashMap<String, AlignmentContextBySample> contexts, StratifiedContext contextType) {
        AlignmentContextBySample context = contexts.get(POOL_SAMPLE_NAME);
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context.getContext(contextType));

        System.out.printf("%s %.2f %d%n", alt, f, pileup.getBases().length());
        // TODO -- IMPLEMENT ME

        return -100.0;
    }
}