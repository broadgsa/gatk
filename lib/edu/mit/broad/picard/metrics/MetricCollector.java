package edu.mit.broad.picard.metrics;

import edu.mit.broad.sam.SAMRecord;

/**
 * Interface for objects that collect metrics about SAMRecords.
 */
public interface MetricCollector<T extends MetricBase> {
    T getMetrics();
    
    /** Called after collector is constructed to populate the metrics object. */
    void setMetrics(T metrics);
    
    /**
     * Called when collection is complete. Implementations can do any calculations
     * that must wait until all records are visited at this time.
     */
    void onComplete();

    /**
     * Visitor method called to have the record considered by the collector.
     */
    void addRecord(SAMRecord record);
}