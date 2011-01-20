package org.broadinstitute.sting.utils.threading;

import java.util.concurrent.locks.ReentrantLock;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: 1/19/11
 * Time: 9:50 AM
 *
 * Simple extension of a ReentrantLock that supports a close method
 */
public class ClosableReentrantLock extends ReentrantLock {
    public void close() {}
}
