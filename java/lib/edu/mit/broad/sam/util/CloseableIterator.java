/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2008 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package edu.mit.broad.sam.util;

import java.util.Iterator;

/**
 * This interface is used by iterators that use releasable resources during iteration.
 * 
 * The consumer of a CloseableIterator should ensure that the close() method is always called,
 * for example by putting such a call in a finally block.  Two conventions should be followed
 * by all implementors of CloseableIterator:
 * 1) The close() method should be idempotent.  Calling close() twice should have no effect.
 * 2) When hasNext() returns false, the iterator implementation should automatically close itself.
 *    The latter makes it somewhat safer for consumers to use the for loop syntax for iteration:
 *    for (Type obj : getCloseableIterator()) { ... }
 * 
 * We do not inherit from java.io.Closeable because IOExceptions are a pain to deal with.
 */
public interface CloseableIterator<T>
    extends Iterator<T> {

    public void close();
}
