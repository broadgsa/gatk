package edu.mit.broad.picard.util;

import edu.mit.broad.picard.util.Histogram.Bin;

import java.util.TreeMap;

/**
 * Class for computing and accessing histogram type data.  Stored internally in
 * a sorted Map so that keys can be iterated in order.
 *
 * @author Tim Fennell
 */
public class Histogram<K extends Comparable> extends TreeMap<K, Bin> {
    private String binLabel   = "BIN";
    private String valueLabel = "VALUE";
    private double count = 0;
    private Double mean;

    /** Constructs a new Histogram with default bin and value labels. */ 
    public Histogram() { }

    /** Constructs a new Histogram with supplied bin and value labels. */
    public Histogram(String binLabel, String valueLabel) {
        this.binLabel = binLabel;
        this.valueLabel = valueLabel;
    }

    /** Represents a bin in the Histogram. */
    public class Bin {
        private final K id;
        private double value = 0;

        /** Constructs a new bin with the given ID. */
        private Bin(K id) { this.id = id; }

        /** Gets the ID of this bin. */
        public K getId() { return id; }

        /** Gets the value in the bin. */
        public double getValue() { return value; }

        /** Returns the String format for the value in the bin. */
        public String toString() { return String.valueOf(this.value); }

        /** Checks the equality of the bin by ID and value. */
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Bin bin = (Bin) o;

            if (Double.compare(bin.value, value) != 0) return false;
            if (!id.equals(bin.id)) return false;

            return true;
        }
        
        public double getIdValue() {
            if (id instanceof Number) {
                return ((Number) id).doubleValue();
            } else {
                throw new UnsupportedOperationException("getIdValue only supported for Histogram<? extends Number>");
            }
        }
    }

    /** Prefill the histogram with the supplied set of bins. */
    public void prefillBins(K... ids) {
        for (K id : ids) {
            put(id, new Bin(id));
        }
    }

    /** Increments the value in the designated bin by 1. */
    public void increment(K id) {
        increment(id, 1d);
    }

    /** Increments the value in the designated bin by the supplied increment. */
    public void increment(K id, double increment) {
        Bin bin = get(id);
        if (bin == null) {
            bin = new Bin(id);
            put(id, bin);
        }

        bin.value += increment;
        count += increment;
        mean = null;
    }

    public String getBinLabel() { return binLabel; }
    public void setBinLabel(String binLabel) { this.binLabel = binLabel; }

    public String getValueLabel() { return valueLabel; }
    public void setValueLabel(String valueLabel) { this.valueLabel = valueLabel; }

    /** Checks that the labels and values in the two histograms are identical. */
    public boolean equals(Object o) {
        return o != null &&
                (o instanceof Histogram) &&
                ((Histogram) o).binLabel.equals(this.binLabel) &&
                ((Histogram) o).valueLabel.equals(this.valueLabel) &&
                super.equals(o);
    }

    public double getMean() {
        if (mean == null) {
            double total = 0;
            for (Bin bin : values()) {
                total += bin.getValue() * bin.getIdValue();
            }
    
            mean = total / count;
        }
        
        return mean;
    }
    
    public double getStandardDeviation() {
        double total = 0;
        for (Bin bin : values()) {
            total += bin.getValue() * bin.getIdValue() * bin.getIdValue();
        }

        return Math.sqrt((total / count) - (getMean() * getMean()));
    }
    
    public double getMedian() {
        double total = 0;
        double halfCount = count / 2;
        for (Bin bin : values()) {
            total += bin.getValue();
            if (total >= halfCount) {
                return bin.getIdValue();
            }
        }
        return 0;
    }

    public double getMin() {
        return firstEntry().getValue().getIdValue();
    }
    
    public double getMax() {
        return lastEntry().getValue().getIdValue();
    }
    
    public double getCount() {
        return count;
    }
}
