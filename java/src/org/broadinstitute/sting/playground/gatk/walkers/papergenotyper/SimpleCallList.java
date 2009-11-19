package org.broadinstitute.sting.playground.gatk.walkers.papergenotyper;

import org.broadinstitute.sting.utils.StingException;

import java.io.*;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collection;

/**
 * Created by IntelliJ IDEA.
 * User: aaronmckenna
 * Date: Nov 19, 2009
 * Time: 12:50:20 AM
 *
 * A simple class, that dumps the records to disk when we've hit a threshold.
 * This class makes the GATKPaperGenotyper much simpler to take in for the reader.
 */
class SimpleCallList extends AbstractList<SimpleCall> {
        private File outputLocation;
        private ArrayList<SimpleCall> list = new ArrayList<SimpleCall>();
        private int WRITE_LIMIT = 100000;
        public SimpleCallList(File writeTo) {
            outputLocation = writeTo;
        }

        public boolean add(SimpleCall call) {
            boolean added = list.add(call);
            writeIfNeeded();
            return added;
        }

        public boolean addAll(Collection<? extends SimpleCall> otherCalls) {
            boolean added = list.addAll(otherCalls);
            writeIfNeeded();
            return added;
        }

        public void writeIfNeeded() {
            synchronized(list) {
                if (list.size() > WRITE_LIMIT) {
                    try {
                        PrintWriter writer = new PrintWriter(new FileWriter(outputLocation, true));
                        for (SimpleCall call : list) writer.println(call.toString());
                        writer.close();
                    } catch (IOException e) {
                        throw new StingException("Unable to write to file " + outputLocation);
                    }
                    list.clear();
                }
            }
        }
        @Override
        public int size() {
            return list.size();
        }

        @Override
        public SimpleCall get(int index) {
            return list.get(index);
        }

    public void close() {
        WRITE_LIMIT = 0;
        writeIfNeeded();
    }
}

