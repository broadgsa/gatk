package edu.mit.broad.sting.utils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.nio.channels.FileChannel;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Mar 2, 2009
 * Time: 2:25:18 PM
 *
 * This class is intended to track the reading of files composed of records of approximately equivalent
 * size.  It can be used to estimate time to completion, complete read, performance of io, etc.
 *
 */
public class FileProgressTracker<T> implements Iterator<T> {
    private static int DEFAULT_HISTORY_SIZE = 1000;

    private int historySize = DEFAULT_HISTORY_SIZE;
    private int samplingFrequency = 1000;
    private File file;
    private FileChannel channel;
    private ArrayList<Long> history;
    private long nNexts = 0;
    private Iterator<T> it = null;
    private long startTime = -1;
    private int historyI = 0;

    public FileProgressTracker( File file, Iterator<T> it, FileChannel channel, int historySize ) {
        this.file = file;
        this.channel = channel;
        this.it = it;
        this.historySize = historySize;
        this.history = new ArrayList<Long>(Collections.nCopies(historySize, 0L));
        startTime = System.currentTimeMillis();
    }

    public FileProgressTracker( File file, Iterator<T> it, FileChannel channel ) {
        this(file, it, channel, DEFAULT_HISTORY_SIZE);
    }

    // -----------------------------------------------------------------
    //
    // iterator support
    //
    // -----------------------------------------------------------------
    public boolean hasNext() { return it.hasNext(); }

    public T next() {
        T x = it.next();
        if ( nNexts % samplingFrequency == 0 ) {
            inc();
            //printStatus();
        }
        nNexts++;
        return x;
    }

    public void remove () {
        it.remove();
    }

    /**
     * Fundamental operation -- must be called ever time a record is read from the file
     *   Enables the system to track the relationship between file byte offsets and record
     *   sizes. 
     */
    public void inc() {
        int i = historyIndex();
        long pos = getPosition();
        history.set(i, pos);
        historyI++;

//        for ( long x : history ) {
//            System.out.printf("%d ", x);
//        }
//        System.out.printf("%n");
//
//        for ( long x : recordSizes() ) {
//            System.out.printf("%d ", x);
//        }
//        System.out.printf("%n");
    }

    public long nRecordsProcessed() {
        return nNexts;
    }

    public double elapsedTimeInSecs() {
        return (System.currentTimeMillis() - startTime) / 1000.0;
    }
    
    public int historyIndex() {
        return historyIndex(historyI);        
    }
    public int historyIndex(long index) {
        return (int)((index + historySize) % historySize);        
    }

    public int averageRecordSize() {
        return Math.round((int)Utils.average(recordSizes(), Math.min(historyI - 1, historySize)) / samplingFrequency);
    }

    public double processingRate() {
        return nRecordsProcessed() / elapsedTimeInSecs();
    }

    public long estRecordsInFile() {
        return (long)(getFileSize() / Math.max(averageRecordSize(),1));
    }

    public double estFractionProgressThroughFile() {
        return (1.0 * nRecordsProcessed()) / estRecordsInFile();
    }

    public double estTimeTotal() {
        return estRecordsInFile() / processingRate();
    }    

    public double estTimeRemaining() {
        return estTimeTotal() * ( 1 - estFractionProgressThroughFile() );
    }

    public void printStatus() {
        System.out.printf("FileProgressTracker:%n");
        System.out.printf("  -> File size is:                                   %d%n", getFileSize());
        System.out.printf("  -> Sampling depth:                                 %d%n", historyI);
        System.out.printf("  -> File position:                                  %d%n", getPosition());
        System.out.printf("  -> Number of records processed:                    %d%n", nRecordsProcessed());
        System.out.printf("  -> Average record size is                          %d%n", averageRecordSize());
        System.out.printf("  -> Elapsed time in secs is                         %.2f%n", elapsedTimeInSecs());
        System.out.printf("  -> Processing rate (records per second)            %.2f%n", processingRate());
        System.out.printf("  -> Estimated number of records in file             %d%n", estRecordsInFile());
        System.out.printf("  -> Estimated percent progress through file         %.2f%n", estFractionProgressThroughFile() * 100.0);
        System.out.printf("  -> Estimated time for entire processing            %.2f hrs / %.2f min / %.2f sec%n", estTimeTotal() / (60*60), estTimeTotal() / (60), estTimeTotal());
        System.out.printf("  -> Estimated time remaining                        %.2f hrs / %.2f min / %.2f sec%n", estTimeRemaining() / (60*60), estTimeRemaining() / 60, estTimeRemaining());
    }

    public String progressMeter() {
        //printStatus();
        return String.format("Est. %.2f%% completed, time remaining (%.2f hrs / %.2f min) of (%.2f hrs / %.2f min) total",
                estFractionProgressThroughFile() * 100.0,
                estTimeRemaining() / (60*60), estTimeRemaining() / 60,
                estTimeTotal() / (60*60), estTimeTotal() / (60));                
    }

    public ArrayList<Long> recordSizes() {
        ArrayList<Long> sizes = new ArrayList<Long>(history);
        for ( int i = 0; i < historySize; i++ ) {
            long val = 0; 
            if ( history.get(historyIndex(i)) >= history.get(historyIndex(i-1) ) )
                val = history.get(historyIndex(i)) - history.get(historyIndex(i-1));
            sizes.set(i, val);
        }

//        for ( long size : sizes ) {
//            System.out.printf("%d ", size);
//        }
//        System.out.printf("%n");

        return sizes;
    }

    private final long getPosition() {
        try {
            return channel.position();
        } catch ( IOException e ) {
            return 0;
        }
    }
    
    private final long getFileSize() {
        return file.length();
    }
    

}
