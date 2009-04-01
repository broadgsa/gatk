package org.broadinstitute.sting.gatk.walkers;

import java.io.PrintStream;

import org.broadinstitute.sting.gatk.GenomeAnalysisTK;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Mar 17, 2009
 * Time: 1:53:31 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class Walker<MapType, ReduceType> {
    // TODO: Can a walker be templatized so that map and reduce live here?

    /**
     * A stream for writing normal (non-error) output.  System.out by default.
     */
    protected PrintStream out = null;

    /**
     * A stream for writing error output.  System.err by default.
     */
    protected PrintStream err = null;

    protected Walker() {
	if( GenomeAnalysisTK.Instance != null ) {
	    GenomeAnalysisTK.Instance.loadArgumentsIntoObject(this);
	    out = GenomeAnalysisTK.Instance.out;
	    err = GenomeAnalysisTK.Instance.err;
	}
	else {
	    out = System.out;
	    err = System.err;
	}
    }

    /**
     * Retrieve the toolkit, for peering into internal structures that can't
     * otherwise be read.  Use sparingly, and discuss uses with software engineering
     * team.
     * @return The genome analysis toolkit.
     */
    protected GenomeAnalysisTK getToolkit() {
        return GenomeAnalysisTK.Instance;
    }

    public void initialize() { }
    public void onTraversalDone(ReduceType result) {
        out.println("[REDUCE RESULT] Traversal result is: " + result);
    }
}
