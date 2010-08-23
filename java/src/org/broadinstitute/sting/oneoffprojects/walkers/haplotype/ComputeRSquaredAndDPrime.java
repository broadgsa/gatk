package org.broadinstitute.sting.oneoffprojects.walkers.haplotype;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;

import java.io.PrintStream;

/**
 * A very simple code snippet from Steve Schaffner that computs R^2 and D' given observed counts for
 * AA, Aa, aA, and aa genotypes.  This code is meant only to be instructive.  Do not use this for
 * anything, ever!
 */
public class ComputeRSquaredAndDPrime extends RodWalker<Integer, Integer> {
    @Argument(fullName="AA", shortName="AA", doc="Number of counts for AA genotype") public Integer AA;
    @Argument(fullName="Aa", shortName="Aa", doc="Number of counts for Aa genotype") public Integer Aa;
    @Argument(fullName="aA", shortName="aA", doc="Number of counts for aA genotype") public Integer aA;
    @Argument(fullName="aa", shortName="aa", doc="Number of counts for aa genotype") public Integer aa;

    @Output
    private PrintStream out;

    public void initialize() {
        int i, j;
        Integer[][] hap = new Integer[2][2];
        for (int k = 0; k < 2; k++) {
            hap[k] = new Integer[2];
        }

        hap[0][0] = AA;
        hap[0][1] = Aa;
        hap[1][0] = aA;
        hap[1][1] = aa;

        Integer[] colTot = new Integer[2];
        Integer[] rowTot = new Integer[2];
        double prod, dprime, f1, f2, ddenom, r2;

        for (i = 0; i < 2; i++) { rowTot[i] = colTot[i] = 0; }
        for (i = 0; i < 2; i++) {
            for (j = 0; j < 2; j++) {
                rowTot[j] += hap[i][j];
                colTot[i] += hap[i][j];
            }
        }
        prod = rowTot[0] * rowTot[1] * colTot[0] * colTot[1];

        if (prod == 0) {
            out.println("Missing data");
            System.exit(1);
        }
        dprime = hap[0][0] * hap[1][1] - hap[0][1] * hap[1][0];
        if (dprime > 0) {
            f1 = rowTot[0] * colTot[1];
            f2 = rowTot[1] * colTot[0];
            ddenom = (f1 > f2) ? f2 : f1;
        } else {
            f1 = rowTot[0] * colTot[0];
            f2 = rowTot[1] * colTot[1];
            ddenom = (f1 > f2) ? f2 : f1;
        }
        r2 = dprime * dprime / prod;
        dprime /= ddenom;

        out.printf("r2: %.5f%n", r2);
        out.printf("D': %.5f%n", dprime);

        System.exit(0);
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return null;
    }

    @Override
    public Integer reduceInit() {
        return null;
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        return null;
    }
}
