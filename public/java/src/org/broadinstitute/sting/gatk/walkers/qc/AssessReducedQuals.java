package org.broadinstitute.sting.gatk.walkers.qc;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.io.PrintStream;

/**
 * Created with IntelliJ IDEA.
 * User: ami
 * Date: 10/19/12
 * Time: 9:09 AM
 * To change this template use File | Settings | File Templates.
 */
public class AssessReducedQuals extends LocusWalker<GenomeLoc, GenomeLoc> implements TreeReducible<GenomeLoc> {

    private static final String original = "original";
    private static final String reduced = "reduced";

    @Output
    protected PrintStream out;

    @Override
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    public void initialize() {}    //todo: why we need that?

    @Override
    public GenomeLoc map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return null;

        int originalQualsIndex = 0;
        int reducedQualsIndex = 1;
        double epsilon = 0;
        double[] quals = getPileupQuals(context.getBasePileup());
        return (quals[originalQualsIndex] - quals[reducedQualsIndex] >= epsilon || quals[originalQualsIndex] - quals[reducedQualsIndex] <= -1*epsilon) ? ref.getLocus() : null;

    }

    private double[] getPileupQuals(final ReadBackedPileup readPileup) {

        int originalQualsIndex = 0;
        int reducedQualsIndex = 1;
        double[] quals = new double[2];

        for( PileupElement p : readPileup){
            if ( (int)p.getQual() > 2 && p.getMappingQual() > 0 && !p.isDeletion() ){
                if (p.getRead().isReducedRead()){
                    double tempQual = (double)(p.getQual()) * p.getRepresentativeCount();
                    quals[reducedQualsIndex] += tempQual;
                }
                else
                {
                    double tempQual = (double)(p.getQual());
                    quals[originalQualsIndex] += tempQual;
                }
            }

        }


        return quals;
    }

    //public void onTraversalDone(GenomeLoc sum) {
    //    if ( sum != null )
    //        out.println(sum);
    //}

    @Override
    public GenomeLoc treeReduce(GenomeLoc lhs, GenomeLoc rhs) {
        if ( lhs == null )
            return rhs;

        if ( rhs == null )
            return lhs;

        // if contiguous, just merge them
        if ( lhs.contiguousP(rhs) )
            return getToolkit().getGenomeLocParser().createGenomeLoc(lhs.getContig(), lhs.getStart(), rhs.getStop());

        // otherwise, print the lhs and start over with the rhs
        out.println(lhs);
        return rhs;
    }

    @Override
    public GenomeLoc reduceInit() {
        return null;
    }

    @Override
    public GenomeLoc reduce(GenomeLoc value, GenomeLoc sum) {
        if ( value == null )
            return sum;

        if ( sum == null )
            return value;

        // if contiguous, just merge them
        if ( sum.contiguousP(value) )
            return getToolkit().getGenomeLocParser().createGenomeLoc(sum.getContig(), sum.getStart(), value.getStop());

        // otherwise, print the sum and start over with the value
        out.println(sum);
        return value;
    }
}
