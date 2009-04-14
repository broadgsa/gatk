package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodSAMPileup;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Utils;
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
    @Argument(fullName="verbose",required=false,defaultValue="false")
    public boolean VERBOSE;

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);

        rodSAMPileup truePileup = (rodSAMPileup)tracker.lookup("pileup", null);
        if ( truePileup == null )
            Utils.scareUser(String.format("No pileup data available at %s given GATK's output of %s -- this walker requires samtools pileup data over all bases",
                    context.getLocation(), pileup.getBases()));

        String pileupDiff = BasicPileup.pileupDiff(pileup, truePileup);
        if ( pileupDiff != null ) {
            out.printf("%s vs. %s%n", pileup.getPileupString(), truePileup.getPileupString());
            throw new RuntimeException(String.format("Pileups aren't equal: %s", pileupDiff));
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