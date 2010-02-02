package org.broadinstitute.sting.gatk.walkers.varianteval;

import java.util.ArrayList;

/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */
public class Histogram<T> {
    ArrayList<T> data;
    int nBins = 100;
    double minValue, maxValue;

    public Histogram(int nBins, double minValue, double maxValue, T initialValue) {
        data = new ArrayList<T>(nBins);
        this.nBins = nBins;
        this.minValue = minValue;
        this.maxValue = maxValue;
        initialize(initialValue);
    }

    public void initialize(T initialValue) {
        data.clear();
        for ( int i = 0; i < nBins; i++ ) {
            data.add(i, initialValue);
        }
    }

    public void setBin(int i, T value) {
        data.set(i, value);
    }

    public void setBin(double x, T value) {
        setBin(x2bin(x), value);
    }

    public T getBin(int binI) {
        return data.get(binI);
    }

    public T getBin(double x) {
        return getBin(x2bin(x));
    }

    public int x2bin(double x) {
        return (int)Math.floor((x - minValue) / ((maxValue - minValue) * nBins));
    }

    public double bin2x(int bin) {
        return (maxValue - minValue) * ((bin + 0.5) / nBins) + minValue;
    }
}
