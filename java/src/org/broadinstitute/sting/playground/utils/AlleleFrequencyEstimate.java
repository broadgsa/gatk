package org.broadinstitute.sting.playground.utils;

import org.broadinstitute.sting.playground.gatk.walkers.AlleleFrequencyWalker;

public class AlleleFrequencyEstimate {
    //AlleleFrequencyEstimate();
    public String location;
    public char ref;
    public char alt;
    public int N;
    public double qhat;
    public double qstar;
    public double LOD;
    public int depth;

    public double getQstar() 
    {
        return qstar;
    }

    public double getLOD() 
    {
        return LOD;
    }

    public AlleleFrequencyEstimate(String location, char ref, char alt, int N, double qhat, double qstar, double LOD, int depth)
    {
        this.location = location;
        this.ref = ref;
        this.alt = alt;
        this.N = N;
        this.qhat = qhat;
        this.qstar = qstar;
        this.LOD = LOD;
        this.depth = depth;
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
