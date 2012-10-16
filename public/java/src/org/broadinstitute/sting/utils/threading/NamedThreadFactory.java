package org.broadinstitute.sting.utils.threading;

import java.util.concurrent.ThreadFactory;

/**
 * Thread factor that produces threads with a given name pattern
 *
 * User: depristo
 * Date: 9/5/12
 * Time: 9:22 PM
 *
 */
public class NamedThreadFactory implements ThreadFactory {
    static int id = 0;
    final String format;

    public NamedThreadFactory(String format) {
        this.format = format;
        String.format(format, id); // test the name
    }

    @Override
    public Thread newThread(Runnable r) {
        return new Thread(r, String.format(format, id++));
    }
}
