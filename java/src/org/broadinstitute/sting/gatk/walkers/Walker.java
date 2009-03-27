package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.GenomeAnalysisTK;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Mar 17, 2009
 * Time: 1:53:31 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class Walker<ReduceType> {
    // TODO: Can a walker be templatized so that map and reduce live here?

    protected Walker() {
        GenomeAnalysisTK.Instance.loadArgumentsIntoObject(this);
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
        // TODO: replace with the correct output stream
        System.out.println("[REDUCE RESULT] Traversal result is: " + result);
    }
}
