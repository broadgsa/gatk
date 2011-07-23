package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.variantcontext.Genotype;

class NewCompEvalGenotypes {
    private GenomeLoc loc;
    private Genotype compGt;
    private Genotype evalGt;

    public NewCompEvalGenotypes(GenomeLoc loc, Genotype compGt, Genotype evalGt) {
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
