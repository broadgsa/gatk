/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/

package edu.mit.broad.picard.metrics;

import edu.mit.broad.sam.SAMRecord;

public class AggregateMetricCollector<T extends MetricBase> implements MetricCollector<T> {
    private final MetricCollector<T>[] collectors;

    public AggregateMetricCollector(MetricCollector<T>... collectors) {
        if (collectors.length == 0) {
            throw new IllegalArgumentException("Must supply at least one collector.");
        }
        this.collectors = collectors;
    }

    @Override
    public void addRecord(SAMRecord record) {
        for (MetricCollector<T> collector : this.collectors) {
            collector.addRecord(record);
        }
    }

    @Override
    public void onComplete() {
        for (MetricCollector<T> collector : this.collectors) {
            collector.onComplete();
        }
    }

    @Override
    public void setMetrics(T metrics) {
        for (MetricCollector<T> collector : this.collectors) {
            collector.setMetrics(metrics);
        }
    }
    
    @Override
    public T getMetrics() {
        return this.collectors[0].getMetrics();
    }
}