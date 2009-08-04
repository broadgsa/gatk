package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.io.FileNotFoundException;
import java.io.PrintWriter;

/**
 * Split up two call sets into their various s"Venn diagram" sets
 */
@Requires(value={DataSource.REFERENCE},referenceMetaData={@RMD(name="callset1",type=AllelicVariant.class),@RMD(name="callset2",type=AllelicVariant.class)})
public class VennCallSetsWalker extends RefWalker<Integer, Integer> {
    @Argument(fullName="union_out", shortName="union", doc="File to which union of sets should be written", required=false)
    String UNION_OUT = null;
    @Argument(fullName="intersect_out", shortName="intersect", doc="File to which intersection of sets should be written", required=false)
    String INTERSECTION_OUT = null;
    @Argument(fullName="set1_unique_out", shortName="set1", doc="File to which calls unique to set1 should be written", required=false)
    String SET1_OUT = null;
    @Argument(fullName="set2_unique_out", shortName="set2", doc="File to which calls unique to set2 should be written", required=false)
    String SET2_OUT = null;

    private PrintWriter union_writer = null;
    private PrintWriter intersect_writer = null;
    private PrintWriter set1_writer = null;
    private PrintWriter set2_writer = null;

    /**
     * Prepare the output file and the list of available features.
     */
    public void initialize() {
        try {
            if ( UNION_OUT != null ) {
                union_writer = new PrintWriter(UNION_OUT);
            }
            if ( INTERSECTION_OUT != null ) {
                intersect_writer = new PrintWriter(INTERSECTION_OUT);
            }
            if ( SET1_OUT != null ) {
                set1_writer = new PrintWriter(SET1_OUT);
            }
            if ( SET2_OUT != null ) {
                set2_writer = new PrintWriter(SET2_OUT);
            }
        } catch (FileNotFoundException e) {
            throw new StingException(String.format("Could not open file(s) for writing"));
        }
    }

    /**
     * Initialize the number of loci processed to zero.
     *
     * @return 0
     */
    public Integer reduceInit() { return 0; }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        AllelicVariant call1 = (AllelicVariant)tracker.lookup("callset1", null);
        AllelicVariant call2 = (AllelicVariant)tracker.lookup("callset2", null);

        // Ignore places where we don't have a variant or where the reference base is ambiguous.
        if ( union_writer != null )
            printVariant(union_writer, call1 != null ? call1 : call2);
        if ( intersect_writer != null )
            printVariant(intersect_writer, call1 != null && call2 != null ? call1 : null);
        if ( set1_writer != null )
            printVariant(set1_writer, call1 != null && call2 == null ? call1 : null);
        if ( set2_writer != null )
            printVariant(set2_writer, call1 == null && call2 != null ? call2 : null);

        return 1;
    }

    static void printVariant(PrintWriter writer, AllelicVariant variant) {
        if ( variant != null )
            writer.println(variant.toString());
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
     * @param result  the number of loci seen.
     */
    public void onTraversalDone(Integer result) {
        out.printf("Processed %d loci.\n", result);

        if ( union_writer != null )
            union_writer.close();
        if ( intersect_writer != null )
            intersect_writer.close();
        if ( set1_writer != null )
            set1_writer.close();
        if ( set2_writer != null )
            set2_writer.close();
    }
}