package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.gatk.walkers.AlleleFrequencyWalker;

public class AlleleFrequencyEstimate {
    //AlleleFrequencyEstimate();
    char ref;
    char alt;
    int N;
    double qhat;
    double qstar;
    double logOddsVarRef;

    public double getQstar() {
        return qstar;
    }

    public double getLogOddsVarRef() {
        return logOddsVarRef;
    }


    public AlleleFrequencyEstimate(char ref, char alt, int N, double qhat, double qstar, double logOddsVarRef) {
        this.ref = ref;
        this.alt = alt;
        this.N = N;
        this.qhat = qhat;
        this.qstar = qstar;
        this.logOddsVarRef = logOddsVarRef;
    }

    public String asString() {
        // Print out the called bases
        // Notes: switched from qhat to qstar because qhat doesn't work at n=1 (1 observed base) where having a single non-ref
        //        base has you calculate qstar = 0.0 and qhat = 0.49 and that leads to a genotype predicition of AG according
        //        to qhat, but AA according to qstar.  This needs to be further investigated to see whether we really want
        //        to use qstar, but make N (number of chormosomes) switch to n (number of reads at locus) for n=1
        long numNonrefBases = Math.round(qstar * N);
        long numRefBases = N - numNonrefBases;
        return AlleleFrequencyWalker.repeat(ref, numRefBases) + AlleleFrequencyWalker.repeat(alt, numNonrefBases);
    }
}