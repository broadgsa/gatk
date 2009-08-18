package org.broadinstitute.sting.playground.gatk.walkers.indels;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.cmdLine.Argument;

/**
 * filter an indel callset based on given criteria
 */
@Requires(value={DataSource.REFERENCE},referenceMetaData={@RMD(name="indels",type=AllelicVariant.class)})
@Reference(window=@Window(start=-20,stop=20))
public class IndelFilterWalker extends RefWalker<Integer, Integer> {
    @Argument(fullName="homopolymerRunMax", shortName="homopolMax", doc="filter indels within homopolymer runs greater than the given length (max 20)", required=false)
    Integer HOMOPOLYMER_MAX = 20;
    @Argument(fullName="homopolymerRunMin", shortName="homopolMin", doc="filter indels within homopolymer runs less than the given length", required=false)
    Integer HOMOPOLYMER_MIN = 0;
    @Argument(fullName="sizeMax", shortName="sizeMax", doc="filter indels greater than a certain size", required=false)
    Integer SIZE_MAX  = 100;
    @Argument(fullName="sizeMin", shortName="sizeMin", doc="filter indels less than a certain size", required=false)
    Integer SIZE_MIN  = 0;

    public void initialize() {
        if  ( HOMOPOLYMER_MAX > 20 )
            HOMOPOLYMER_MAX = 20;
    }

    public Integer reduceInit() { return 0; }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        AllelicVariant indel = (AllelicVariant)tracker.lookup("indels", null);

        if ( indel == null || !indel.isIndel() )
            return 0;

        if ( indel.length() < SIZE_MIN || indel.length() > SIZE_MAX )
            return 0;
        
        int homopol = homopolymerRunSize(ref, indel);
        if ( homopol < HOMOPOLYMER_MIN || homopol > HOMOPOLYMER_MAX )
            return 0;

        out.println(indel);
        return 1;
    }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {
        out.printf("output %d indels.\n", result);
    }

    private int homopolymerRunSize(ReferenceContext ref, AllelicVariant indel) {
        char[] bases = ref.getBases();
        GenomeLoc window = ref.getWindow();
        GenomeLoc locus = ref.getLocus();

        int refBasePos = (int)(locus.getStart() - window.getStart());
        char indelBase = indel.isDeletion() ? bases[refBasePos+1] : indel.getAltBasesFWD().charAt(0);
        int leftRun = 0;
        for ( int i = refBasePos; i >= 0; i--) {
            if ( bases[i] != indelBase )
                break;
            leftRun++;
        }

        indelBase = indel.isDeletion() ? bases[refBasePos+indel.length()] : indel.getAltBasesFWD().charAt(indel.getAltBasesFWD().length()-1);
        int rightRun = 0;
        for ( int i = refBasePos + (indel.isDeletion() ? 1+indel.length() : 1); i < bases.length; i++) {
            if ( bases[i] != indelBase )
                break;
            rightRun++;
        }

        System.out.println(String.valueOf(bases) + ": " + leftRun + " / " + rightRun);
        return Math.max(leftRun, rightRun);
    }
}