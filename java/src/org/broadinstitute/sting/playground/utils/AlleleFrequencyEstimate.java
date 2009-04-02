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
    public double lodVsRef;
    public double lodVsNextBest;
    public int depth;
    public String notes;

    public AlleleFrequencyEstimate(String location, char ref, char alt, int N, double qhat, double qstar, double lodVsRef, double lodVsNextBest, int depth)
    {
        this.location = location;
        this.ref = ref;
        this.alt = alt;
        this.N = N;
        this.qhat = qhat;
        this.qstar = qstar;
        this.lodVsRef = lodVsRef;
        this.lodVsNextBest = lodVsNextBest;
        this.depth = depth;
        this.notes = "";
    }

    public String asGFFString()
    {
        String[] tokens;
        tokens = location.split(":");
        return String.format("%s\tCALLER\tVARIANT\t%s\t%s\t%f\t.\t.\tREF %c\t;\tALT %c\t;\tFREQ %f\n",
                               tokens[0],
                               tokens[1],
                               tokens[1],
                               lodVsRef,
                               ref,
                               alt,
                               qhat);
    }

    public String asTabularString() {
        return String.format("RESULT %s %c %c %f %f %f %f %d %s\n",
	                                        location,
	                                        ref,
	                                        alt,
	                                        qhat,
	                                        qstar,
                                            lodVsRef,
                                            lodVsNextBest,
	                                        depth, 
                                            notes);
    }

    public String toString() { return asTabularString(); }

    public String asString() {
        // Print out the called bases
        // Notes: switched from qhat to qstar because qhat doesn't work at n=1 (1 observed base) where having a single non-ref
        //        base has you calculate qstar = 0.0 and qhat = 0.49 and that leads to a genotype predicition of AG according
        //        to qhat, but AA according to qstar.  This needs to be further investigated to see whether we really want
        //        to use qstar, but make N (number of chormosomes) switch to n (number of reads at locus) for n=1
        long numNonrefBases = Math.round(qstar * N);
        long numRefBases = N - numNonrefBases;
        if (ref < alt) { // order bases alphabetically
            return AlleleFrequencyWalker.repeat(ref, numRefBases) + AlleleFrequencyWalker.repeat(alt, numNonrefBases);
        }else{
            return AlleleFrequencyWalker.repeat(alt, numNonrefBases) + AlleleFrequencyWalker.repeat(ref, numRefBases);
        }
    }
}
