package org.broadinstitute.sting.gatk.walkers;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Mar 17, 2009
 * Time: 1:53:31 PM
 * To change this template use File | Settings | File Templates.
 */
public interface Walker {
    // TODO: Can a walker be templatized so that map and reduce live here?
    String getName();
    void initialize();
    void onTraversalDone();
}
