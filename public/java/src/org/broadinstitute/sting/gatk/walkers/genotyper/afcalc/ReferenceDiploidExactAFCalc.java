package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

public class ReferenceDiploidExactAFCalc extends DiploidExactAFCalc {
    protected ReferenceDiploidExactAFCalc(int nSamples, int maxAltAlleles, int maxAltAllelesForIndels, final int ploidy) {
        super(nSamples, maxAltAlleles, maxAltAllelesForIndels, ploidy);
    }
}
