package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

public class ReferenceDiploidExactAFCalc extends DiploidExactAFCalc {
    protected ReferenceDiploidExactAFCalc(int nSamples, int maxAltAlleles, final int ploidy) {
        super(nSamples, maxAltAlleles, ploidy);
    }
}
