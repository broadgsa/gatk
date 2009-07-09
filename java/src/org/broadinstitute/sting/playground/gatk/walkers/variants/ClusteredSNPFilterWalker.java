package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodVariants;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.playground.utils.AlleleFrequencyEstimate;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.io.*;
import java.util.*;

/**
 * ClusteredSNPFilterWalker takes a list of variant sites and filters out those that are
 * too clustered together.  At the moment, the variants are expected to be in gelitext format.
 */
@Requires(value={DataSource.REFERENCE},referenceMetaData=@RMD(name="variant",type=rodVariants.class))
public class ClusteredSNPFilterWalker extends RefWalker<Integer, Integer> {
    @Argument(fullName="windowSize", shortName="window", doc="window size for calculating clusters", required=false)
    int windowSize = 10;
    @Argument(fullName="clusterSize", shortName="cluster", doc="number of variants in a window to be considered a cluster", required=false)
    int clusterSize = 3;
    @Argument(fullName="variants_out", shortName="VO", doc="File to which variants passing the filter should be written", required=true)
    File VARIANTS_OUT = null;
    @Argument(fullName="failed_out", shortName="FO", doc="File to which variants failing the filter should be written", required=false)
    File FILTERED_OUT = null;

    private PrintWriter vwriter = null;
    private PrintWriter fwriter = null;

    private LinkedList<VariantState> queue = new LinkedList<VariantState>();

    /**
     * Prepare the output file and the list of available features.
     */
    public void initialize() {
        try {
            vwriter = new PrintWriter(VARIANTS_OUT);
            vwriter.println(AlleleFrequencyEstimate.geliHeaderString());
            if ( FILTERED_OUT != null ) {
                fwriter = new PrintWriter(FILTERED_OUT);
                fwriter.println(AlleleFrequencyEstimate.geliHeaderString());
            }
        } catch (FileNotFoundException e) {
            throw new StingException(String.format("Could not open file for writing"));
        }
    }

    /**
     * Initialize the number of loci processed to zero.
     *
     * @return 0
     */
    public Integer reduceInit() { return 0; }

    /**
     * For each site of interest, rescore the genotype likelihoods by applying the specified feature set.
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return 1 if the locus was successfully processed, 0 if otherwise
     */
    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        rodVariants variant = (rodVariants)tracker.lookup("variant", null);

        if (variant != null ) {
            // add the new variant to the queue
            queue.offer(new VariantState(variant, context.getLocation()));

            if ( queue.size() < clusterSize )
                return 1;

            // remove the head of the queue if applicable
            if ( queue.size() > clusterSize ) {
                VariantState var = queue.remove();
                if ( var.passed )
                    vwriter.println(var.variant);
                else if ( fwriter != null)
                    fwriter.println(var.variant);
            }

            VariantState head = queue.peek();
            if ( context.getLocation().getContigIndex() == head.location.getContigIndex() &&
                 Math.abs(head.location.getStart() - context.getLocation().getStart()) <= windowSize ) {
                for (int i = 0; i < clusterSize; i++)
                    queue.get(i).passed = false;
            }
        }

        return 1;
    }

    /**
     * Increment the number of loci processed.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return the new number of loci processed.
     */
    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    /**
     * Tell the user the number of loci processed and close out the new variants file.
     *
     * @param result  the number of variants seen.
     */
    public void onTraversalDone(Integer result) {
        while ( queue.size() > 0 ) {
            VariantState var = queue.remove();
            if ( var.passed )
                vwriter.println(var.variant);
            else if ( fwriter != null)
                fwriter.println(var.variant);
        }

        out.printf("Processed %d variants.\n", result);

        vwriter.close();
        if ( fwriter != null )
            fwriter.close();
    }

    private class VariantState {
        public rodVariants variant;
        public GenomeLoc location;
        public boolean passed = true;

        public VariantState(rodVariants variant, GenomeLoc location) {
            this.variant = variant;
            this.location = location;
        }
    }
}