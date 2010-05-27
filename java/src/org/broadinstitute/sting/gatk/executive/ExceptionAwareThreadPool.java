package org.broadinstitute.sting.gatk.executive;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;

import java.util.concurrent.*;

/**
 * an override of the ThreadedPoolExecutor, that throws the exception when an exception is seen in a thread.
 */
public class ExceptionAwareThreadPool extends ThreadPoolExecutor {
    /**
     * our log, which we want to capture anything from this class
     */
    private static Logger logger = Logger.getLogger(ExceptionAwareThreadPool.class);

    public ExceptionAwareThreadPool(int numberOfThreads) {
        super(numberOfThreads, numberOfThreads, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());
        Utils.warnUser("Using the etd (enable threaded debugging) mode and the ExceptionAwareThreadPool is dangerous to " +
                       "yourself and others. DONT USE IT IN A PRODUCTION SETTING, OR AT ALL!!!!");
    }

    /**
     * attempt to determine the fate of a runnable object
     * @param r the runnable (in our case a Future object)
     * @param t any throwables from the thread.
     */
    @Override
    public void afterExecute(Runnable r, Throwable t) {
        super.afterExecute(r, t);
        // from http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=6459119.  The throwable in the method parameters
        // is not actually the thrown exception from the underlying thread, we have to go get the cause.
        if (r instanceof Future<?>) {
            try {
                Object result = ((Future<?>) r).get();
                // once we have the result, we know everything went fine
                logger.debug("Thread completed successfully");
            } catch (InterruptedException ie) {
                caughtException(ie);
            } catch (ExecutionException ee) {
                caughtException(ee.getCause());
            } catch (CancellationException ce) {
                caughtException(ce);
            }
        }
    }

    /**
     * ungracefully shutdown the GATK
     * @param e the throwable object, which caused us to fail
     */
    public void caughtException(Throwable e) {
        // shutdown all the threads, not waiting to finish
        this.shutdownNow();

        // cite the reason we crashed out
        logger.fatal("Thread pool caught an exception from a thread: ");
        e.printStackTrace();

        // bail in the ugliest way possible
        System.exit(1);
    }

}
