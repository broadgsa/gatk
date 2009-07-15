package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.refdata.AllelicVariant;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.LocusContext;

import java.io.PrintStream;
import java.util.List;
import java.util.ArrayList;

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
    protected PrintStream out;
    protected VariantEvalWalker master;
    protected String filename;

    public BasicVariantAnalysis(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public List<String> getParams() {
        return new ArrayList<String>();
    }

    public void initialize(VariantEvalWalker master, PrintStream out, String filename) {
        this.master = master;
        this.out = out;
        this.filename = filename;
    }

    public PrintStream getSummaryPrintStream() {
        return out;
    }

    public PrintStream getCallPrintStream() {
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