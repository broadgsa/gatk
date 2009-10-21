package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.io.PrintStream;
import java.util.List;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Sep 3, 2009
 * Time: 5:04:40 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class BasicPoolVariantAnalysis implements VariantAnalysis{
    protected int numIndividualsInPool;
    protected String name;
    protected String filenames;
    protected PrintStream out, callOut;
    private VariantEvalWalker master;

    public BasicPoolVariantAnalysis(String name, int nIndividuals) {
        this.name = name;
        this.numIndividualsInPool = nIndividuals;
    }

    public String getName() { return this.name; }

    public int getNumberOfIndividualsInPool() { return this.numIndividualsInPool; }

    public PrintStream getSummaryPrintStream() { return this.out; }

    public PrintStream getCallPrintStream() { return this.callOut; }

    public VariantEvalWalker getMaster() { return this.master; }

    public List<String> getParams() { return new ArrayList<String>(); }

    public List<String> done() { return new ArrayList<String>(); }

    public int getNumberOfAllelesInPool() { return 2*getNumberOfIndividualsInPool(); }

    public void initialize(VariantEvalWalker master, PrintStream out, PrintStream callOut, String filenames) {
        this.master = master;
        this.out = out;
        this.callOut = callOut;
        this.filenames = filenames;
    }

    public void finalize(long nSites) {
        // no need to finalize data in general
    }

    public abstract String update(Variation variant, RefMetaDataTracker tracker, char ref, AlignmentContext context);

}
