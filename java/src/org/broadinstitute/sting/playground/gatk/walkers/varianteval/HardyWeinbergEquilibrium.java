package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.refdata.AllelicVariant;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.SNPCallFromGenotypes;
import org.broadinstitute.sting.gatk.refdata.PooledEMSNPROD;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.genotype.HardyWeinbergCalculation;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import cern.jet.math.Arithmetic;

/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */
public class HardyWeinbergEquilibrium extends ViolationVariantAnalysis implements PopulationAnalysis {
    private double threshold;
    int nSites = 0;
    int nViolations = 0;

    HardyWeinbergEquilibrium(double threshold) {
        super("hwe");
        this.threshold = threshold;
    }

    public String update(AllelicVariant eval, RefMetaDataTracker tracker, char ref, LocusContext context) {
        String r = null;

        if ( eval != null &&
             eval instanceof SNPCallFromGenotypes ) {
            nSites++;
            SNPCallFromGenotypes call = (SNPCallFromGenotypes)eval;
            //double pPoint = pointHWcalculate(call);
            double pTailed = tailedHWcalculate(call);
            double p = pTailed;

            //System.out.printf("HWE point=%.4e %s vs. tailed=%.4e %s for %s%n", pPoint, pPoint < threshold ? "***" : "   ", pTailed, pTailed < threshold ? "***" : "   ", call);

            if ( p < threshold ) {
                r = String.format("HWE-violation %f < %f at %s", p, threshold, eval);
                nViolations++;
            }
        }

        return r;
    }

    public double tailedHWcalculate(SNPCallFromGenotypes call) {
        int obsAA = call.nHomRefGenotypes();
        int obsAB = call.nHetGenotypes();
        int obsBB = call.nHomVarGenotypes();
        return HardyWeinbergCalculation.hwCalculate(obsAA, obsAB, obsBB);
    }

    public double pointHWcalculate(SNPCallFromGenotypes call) {
        int nAA = call.nHomRefGenotypes();
        int nAa = call.nHetGenotypes();
        int naa = call.nHomVarGenotypes();
        int nA  = 2 * nAA + nAa;
        int n   = nAA + nAa + naa;

        //
        // from Emigh 1980
        //
        // P = pr[nAa | nA] = multinomial[n over nAA, nAa, naa] / binomial[2n over nA] * 2^nAa
        //
        // where nAA, nAa, naa are the observed numbers of the three genotypes, AA, Aa, and aa,
        // respectively, and nA is the number of A alleles, where nA = 2nAA + nAa, and n is the number of alleles
        //
        int[] mXs = { nAA, nAa, naa };                      // counts of each genotype as vector
        double m = MathUtils.multinomial(mXs);
        double b = Arithmetic.binomial(2 * n, nA);
        double tosses = Math.pow(2, nAa);
        double p = (m / b) * tosses;

        if ( false ) {
            System.out.printf("HWE-violation at %s %f < %f %1.2f %5d %5d %5d %5d %5d %.2e %.2e %.2e => %.6e [%s]%n",
                    call.getLocation(), p, threshold, call.getMAF(), nAA, nAa, naa, nA, n, m, b, tosses, p, call);
            System.out.printf("(factorial(%d) / (factorial(%d) * factorial(%d) * factorial(%d))) / choose(%d, %d) * 2^%d - %f < 1e-3%n",
                    nAA + nAa + naa, nAA, nAa, naa, 2 * n, nA, nAa, p);
        }

        return p;
    }

    public List<String> done() {
        List<String> s = new ArrayList<String>();
        s.add(String.format("n_calls              %d", nSites));
        s.add(String.format("n_violations         %d", nViolations));
        s.add(String.format("violations_rate      %.2f", (100.0*nViolations) / nSites));
        return s;
    }
}