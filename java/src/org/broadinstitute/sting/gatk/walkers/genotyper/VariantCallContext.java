package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;

/**
 * Created by IntelliJ IDEA.
 * User: depristo, ebanks
 * Date: Jan 22, 2010
 * Time: 2:25:19 PM
 *
 * Useful helper class to communicate the results of calculateGenotype to framework
 */
public class VariantCallContext {
    public VariantContext vc = null;
    public byte refBase;

    // Was the site called confidently, either reference or variant?
    public boolean confidentlyCalled = false;

    VariantCallContext(VariantContext vc, boolean confidentlyCalledP) {
        this.vc = vc;
        this.confidentlyCalled = confidentlyCalledP;
    }

    VariantCallContext(VariantContext vc, byte ref, boolean confidentlyCalledP) {
        this.vc = vc;
        this.refBase = ref;
        this.confidentlyCalled = confidentlyCalledP;
    }

    // blank variant context => we're a ref site
    VariantCallContext(boolean confidentlyCalledP) {
        this.confidentlyCalled = confidentlyCalledP;
    }

    public void setRefBase(byte ref) {
        this.refBase = ref;
    }
}