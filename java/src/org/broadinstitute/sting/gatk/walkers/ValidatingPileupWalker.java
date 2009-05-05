package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodSAMPileup;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.Pileup;
import org.broadinstitute.sting.utils.BasicPileup;
import org.broadinstitute.sting.utils.ReadBackedPileup;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 3:22:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class ValidatingPileupWalker extends LocusWalker<Integer, ValidationStats> {
    @Argument(fullName="continue_after_error",doc="Continue after an error",required=false,defaultValue="false")
    public boolean CONTINUE_AFTER_AN_ERROR;

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        Pileup pileup = new ReadBackedPileup(ref, context);
        Pileup truePileup = (Pileup)tracker.lookup("pileup", null);
        
        if ( truePileup == null ) {
            System.out.printf("No truth pileup data available at %s%n", pileup.getPileupString());
            if ( ! CONTINUE_AFTER_AN_ERROR ) {
                Utils.scareUser(String.format("No pileup data available at %s given GATK's output of %s -- this walker requires samtools pileup data over all bases",
                                context.getLocation(), pileup.getBases()));
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