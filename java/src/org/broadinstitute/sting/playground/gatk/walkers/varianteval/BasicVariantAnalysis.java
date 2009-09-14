package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */
public abstract class BasicVariantAnalysis implements VariantAnalysis {
    protected String name;
    protected PrintStream out, callOut;
    private VariantEvalWalker master;
    protected String filename;

    public BasicVariantAnalysis(String name) {
        this.name = name;
    }

    public VariantEvalWalker getMaster() { return master; }

    public String getName() {
        return name;
    }

    public List<String> getParams() {
        return new ArrayList<String>();
    }

    public void initialize(VariantEvalWalker master, PrintStream out, PrintStream callOut, String filename) {
        this.master = master;
        this.out = out;
        this.callOut = callOut;
        this.filename = filename;
    }

    public PrintStream getSummaryPrintStream() {
        return out;
    }

    public PrintStream getCallPrintStream() {
        return callOut;
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

    public abstract String update(Variation eval, RefMetaDataTracker tracker, char ref, AlignmentContext context);
}