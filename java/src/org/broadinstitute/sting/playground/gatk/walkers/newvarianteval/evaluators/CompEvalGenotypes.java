package org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.evaluators;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broadinstitute.sting.utils.GenomeLoc;

class CompEvalGenotypes {
    private GenomeLoc loc;
    private Genotype compGt;
    private Genotype evalGt;

    public CompEvalGenotypes(GenomeLoc loc, Genotype compGt, Genotype evalGt) {
        this.loc = loc;
        this.compGt = compGt;
        this.evalGt = evalGt;
    }

    public GenomeLoc getLocus() {
        return loc;
    }

    public Genotype getCompGenotpye() {
        return compGt;
    }
    public Genotype getEvalGenotype() {
        return evalGt;
    }

    public void setCompGenotype(Genotype compGt) {
        this.compGt = compGt;
    }

    public void setEvalGenotype(Genotype evalGt) {
        this.evalGt = evalGt;
    }
}
