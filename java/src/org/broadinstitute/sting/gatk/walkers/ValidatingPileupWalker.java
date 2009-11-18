package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodSAMPileup;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.Pileup;
import org.broadinstitute.sting.utils.BasicPileup;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.StingException;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 3:22:14 PM
 * To change this template use File | Settings | File Templates.
 */
@Requires(value={DataSource.READS,DataSource.REFERENCE},referenceMetaData=@RMD(name="pileup",type=rodSAMPileup.class))
public class ValidatingPileupWalker extends LocusWalker<Integer, ValidationStats>  implements TreeReducible<ValidationStats> {
    @Argument(fullName="continue_after_error",doc="Continue after an error",required=false)
    public boolean CONTINUE_AFTER_AN_ERROR = false;

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        Pileup pileup = new ReadBackedPileup(ref.getBase(), context);
        Pileup truePileup = getTruePileup( tracker );

        if ( truePileup == null ) {
            out.printf("No truth pileup data available at %s%n", pileup.getPileupString());
            if ( ! CONTINUE_AFTER_AN_ERROR ) {
                Utils.scareUser(String.format("No pileup data available at %s given GATK's output of %s -- this walker requires samtools pileup data over all bases",
                                context.getLocation(), pileup.getBasesAsString()));
            }
        } else {
            String pileupDiff = BasicPileup.pileupDiff(pileup, truePileup);
            if ( pileupDiff != null ) {
                out.printf("%s vs. %s%n", pileup.getPileupString(), truePileup.getPileupString());
                if ( ! CONTINUE_AFTER_AN_ERROR ) {
                    throw new RuntimeException(String.format("Pileups aren't equal: %s", pileupDiff));
                }
            }
        }

        return pileup.size();
    }

    // Given result of map function
    public ValidationStats reduceInit() { return new ValidationStats(); }
    public ValidationStats reduce(Integer value, ValidationStats sum) {
        sum.nLoci++;
        sum.nBases += value;
        return sum;
    }

    public ValidationStats treeReduce( ValidationStats lhs, ValidationStats rhs ) {
        ValidationStats combined = new ValidationStats();
        combined.nLoci = lhs.nLoci + rhs.nLoci;
        combined.nBases = lhs.nBases + rhs.nBases;
        return combined;
    }

    /**
     * Extracts the true pileup data from the given rodSAMPileup.  Note that this implementation
     * assumes that the genotype will only be point or indel.
     * @param tracker ROD tracker from which to extract pileup data.
     * @return True pileup data.
     */
    private Pileup getTruePileup( RefMetaDataTracker tracker ) {
        rodSAMPileup pileup = (rodSAMPileup)tracker.lookup("pileup", null);

        if( pileup == null )
            return null;

        if( pileup.hasPointGenotype() )
            return (Pileup)pileup.getPointGenotype();
        else if( pileup.hasIndelGenotype() )
            return (Pileup)pileup.getIndelGenotype();
        else
            throw new StingException("Unsupported pileup type: " + pileup);
    }
}

class ValidationStats {
    public long nLoci = 0;
    public long nBases = 0;

    public ValidationStats() {
    }

    public String toString() {
        return String.format("Validated %d sites covered by %d bases%n", nLoci, nBases);
    }
}