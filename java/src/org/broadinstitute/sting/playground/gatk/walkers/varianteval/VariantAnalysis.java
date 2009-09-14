package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.io.PrintStream;
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
public interface VariantAnalysis {
    public String getName();
    public PrintStream getSummaryPrintStream();
    public PrintStream getCallPrintStream();
    public List<String> getParams();
    public void initialize(VariantEvalWalker master, PrintStream out, PrintStream callOut, String filename);
    public String update(Variation eval, RefMetaDataTracker tracker, char ref, AlignmentContext context);
    public void finalize(long nSites);
    public List<String> done();
}
