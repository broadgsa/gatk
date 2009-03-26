package org.broadinstitute.sting.gatk.walkers;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Mar 17, 2009
 * Time: 1:53:31 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class Walker {
    // TODO: Can a walker be templatized so that map and reduce live here?
    public String getName() {
        // Return name of class, trimming 'Walker' from the end if present.
        String className = getClass().getSimpleName();
        if(className.endsWith(Walker.class.getSimpleName()))
            return className.substring(0,className.lastIndexOf(Walker.class.getSimpleName()));
        else
            return className;
    }

    public void initialize() { }
    public void onTraversalDone() { }
}
