package org.broadinstitute.sting.utils;

import java.io.PrintStream;

/**
 * A useful simple system for timing code.
 *
 * User: depristo
 * Date: Dec 10, 2010
 * Time: 9:07:44 AM
 */
public class SimpleTimer {
    private String name = "";
    private long elapsed = 0l;
    private long startTime = -1l;
    boolean running = false;

    public SimpleTimer(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public SimpleTimer start() {
        elapsed = 0l;
        restart();
        return this;
    }

    public SimpleTimer restart() {
        running = true;
        startTime = currentTime();
        return this;
    }

    public long currentTime() {
        return System.currentTimeMillis();
    }

    public SimpleTimer stop() {
        running = false;
        elapsed += currentTime() - startTime;
        return this;
    }

    public double getElapsedTime() {
        if ( running )
            return (currentTime() - startTime) / 1000.0;
        else
            return elapsed;
    }


    public void printElapsedTime(PrintStream out) {
        out.printf("SimpleTimer %s: %.2f%n", name, getElapsedTime());
    }
}
