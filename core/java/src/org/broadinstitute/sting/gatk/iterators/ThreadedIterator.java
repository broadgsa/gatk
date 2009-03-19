package org.broadinstitute.sting.gatk.iterators;

import java.util.Iterator;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Feb 24, 2009
 * Time: 10:24:38 AM
 * To change this template use File | Settings | File Templates.
 */
public class ThreadedIterator<T> implements Iterator<T>, Runnable {
    private Iterator<T> it;
    private final BlockingQueue<T> queue;
    private int nOps = 0;
    private final int printStateFreq = -1;

    public void run() {
        try {
            while (it.hasNext()) {
                queue.put(it.next());
                printState("addNext");
            }
        } catch (InterruptedException ex) {
            // bail out
            ;
        }
    }

    public synchronized void printState(final String op) {
        if ( printStateFreq != -1 && nOps++ % printStateFreq == 0 )
            System.out.printf("  [%s] Queue has %d elements %d ops%n", op, queue.size(), nOps);
    }

    public ThreadedIterator(Iterator<T> it, int buffSize) {
        this.it = it;
        queue = new LinkedBlockingQueue<T>(buffSize);
        new Thread(this).start();
    }

    public boolean hasNext() {
        return queue.peek() != null || it.hasNext();
    }

    public T next() {
        printState("getNext");
        try {
            return queue.poll(10, TimeUnit.SECONDS);
        } catch (InterruptedException ex) {
            // bail out
            System.out.printf("ThreadedIterator next() timed out...%n");
            printState("getNext");
            return null;
        }
    }

    public void remove () {
        throw new UnsupportedOperationException();
    }
}