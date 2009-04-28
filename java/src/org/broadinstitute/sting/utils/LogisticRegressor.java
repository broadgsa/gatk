package org.broadinstitute.sting.utils;

import net.sf.samtools.SAMFileReader;

import java.util.HashMap;
import java.util.ArrayList;
import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Apr 23, 2009
 * Time: 3:46:30 PM
 * To change this template use File | Settings | File Templates.
 */
public class LogisticRegressor {
    double[][] coefficients;
    int nFeatures;
    int order;
    
    public LogisticRegressor(int nFeatures, int order) {
        this.nFeatures = nFeatures;
        this.order = order;

        if ( nFeatures != 2 )
            throw new IllegalArgumentException("LogisticRegressor currently only supports 2 features :-(");

        // setup coefficient matrix
        coefficients = new double[order+1][order+1];
        for ( int i = 0; i <= order; i++ ) {
            for ( int j = 0; j <= order; j++ ) {
                coefficients[i][j] = 0.0;
            }
        }
    }

    public double[][] getCoefficients() {
        return coefficients;
    }

    public void setCoefficient(int i, int j, double c) {
        coefficients[i][j] = c;
    }

    public double regress(double f1, double f2) {
        double v = 0.0;
        for ( int i = 0; i <= order; i++ ) {
            for ( int j = 0; j <= order; j++ ) {
                double c = coefficients[i][j];
                v += c * Math.pow(f1,i) * Math.pow(f2, j);
                //System.out.printf("i=%d, j=%d, v=%f, c=%f, f1=%f, f2=%f, f1^i=%f, f2^j=%f%n", i, j, v, c, f1, f2, Math.pow(f1,i), Math.pow(f2,j));
            }
        }
        return v;
    }

    public String toString() {
        StringBuilder s = new StringBuilder();
        s.append(String.format("nFeatures=%d, order=%d: ", nFeatures, order));
        for ( int i = 0; i <= order; i++ ) {
            for ( int j = 0; j <= order; j++ ) {
                s.append(" " + coefficients[i][j]);
            }
        }
        
        return s.toString();
    }
}
