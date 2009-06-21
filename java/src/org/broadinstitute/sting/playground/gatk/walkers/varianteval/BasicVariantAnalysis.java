package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.refdata.AllelicVariant;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.LocusContext;

import java.io.PrintStream;
import java.util.List;
import java.util.ArrayList;

public abstract class BasicVariantAnalysis implements VariantAnalysis {
    protected String name;
    protected PrintStream out;
    protected VariantEvalWalker master;

    public BasicVariantAnalysis(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public List<String> getParams() {
        return new ArrayList<String>();
    }

    public void initialize(VariantEvalWalker master, PrintStream out) {
        this.master = master;
        this.out = out;
    }

    public PrintStream getPrintStream() {
        return out;
    }

    public List<String> done() {
        return new ArrayList<String>();
    }

    /**
     * No need to finalize the data in general
     * @param nSites
     */
    public void finalize(long nSites) {

    }

    public abstract String update(AllelicVariant eval, RefMetaDataTracker tracker, char ref, LocusContext context);
}