package org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.evaluators;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.HashMap;

/**
 * Created by IntelliJ IDEA. User: kiran Date: Nov 29, 2010 Time: 3:25:59 PM To change this template use File | Settings
 * | File Templates.
 */
class SamplePreviousGenotypes {
    private HashMap<String, CompEvalGenotypes> sampleGenotypes = null;

    public SamplePreviousGenotypes() {
        this.sampleGenotypes = new HashMap<String, CompEvalGenotypes>();
    }

    public CompEvalGenotypes get(String sample) {
        return sampleGenotypes.get(sample);
    }

    public void put(String sample, CompEvalGenotypes compEvalGts) {
        sampleGenotypes.put(sample, compEvalGts);
    }

    public void put(String sample, GenomeLoc locus, Genotype compGt, Genotype evalGt) {
        sampleGenotypes.put(sample, new CompEvalGenotypes(locus, compGt, evalGt));
    }
}
