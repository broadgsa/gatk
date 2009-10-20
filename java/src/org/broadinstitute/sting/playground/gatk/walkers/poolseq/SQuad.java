package org.broadinstitute.sting.playground.gatk.walkers.poolseq;

/**
 * Created by IntelliJ IDEA.
 * User: Ghost
 * Date: Sep 12, 2009
 * Time: 1:04:40 PM
 * To change this template use File | Settings | File Templates.
 */
public class SQuad<X> extends Quad<X,X,X,X> {
    /* SQuad - single object type Quad
     * or "Simple Quad". Makes it less code-intensive
     * for user to use quads to hold objects.
     */

    public SQuad(X a, X b, X c, X d) { super(a,b,c,d); }
    // don't know why this method even neds to be written
    // but IntelliJ wants it to be there

}
