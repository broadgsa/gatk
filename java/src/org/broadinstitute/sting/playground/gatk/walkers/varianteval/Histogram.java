package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: May 27, 2009
 * Time: 2:37:56 PM
 * To change this template use File | Settings | File Templates.
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
