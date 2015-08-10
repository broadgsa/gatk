/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.tools.walkers.variantrecalibration;

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.text.XReadLines;
import htsjdk.variant.variantcontext.VariantContext;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: 3/11/11
 * Time: 10:33 AM
 * To change this template use File | Settings | File Templates.
 */
public class VQSRCalibrationCurve {
    private final static boolean DEBUG = false;
    List<VQSRRange> points;
    public static final double CERTAIN_FALSE_POSITIVE = -1;

    private static class VQSRRange {
        double start, stop, truePositiveRate;

        public double getStart() {
            return start;
        }

        public double getStop() {
            return stop;
        }

        public double getTruePositiveRate() {
            return truePositiveRate;
        }

        private VQSRRange(double start, double stop, double truePositiveRate) {
            this.start = start;
            this.stop = stop;
            this.truePositiveRate = truePositiveRate;
        }
    }

    public static VQSRCalibrationCurve readFromFile(File source) {
        List<VQSRRange> points = new ArrayList<VQSRRange>();

        try {
            for ( String line : new XReadLines(source).readLines() ) {
                if ( ! line.trim().isEmpty() ) {
                    String[] parts = line.split("\\s+");
                    double fpRate = Double.parseDouble(parts[2]);
                    double tpRate = fpRate >= 1.0 ? CERTAIN_FALSE_POSITIVE : 1.0 - fpRate;
                    points.add(new VQSRRange(Double.parseDouble(parts[0]), Double.parseDouble(parts[1]), tpRate));
                }
            }
        } catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(source, e);
        }

        // ensure that the entire range gets caught
        points.get(0).start = Double.POSITIVE_INFINITY;
        points.get(points.size()-1).stop = Double.NEGATIVE_INFINITY;

        return new VQSRCalibrationCurve(points);
    }

    protected VQSRCalibrationCurve(List<VQSRRange> points) {
        this.points = points;
    }

    public boolean certainFalsePositive(String VQSRQualKey, VariantContext vc) {
        return probTrueVariant(VQSRQualKey, vc) == CERTAIN_FALSE_POSITIVE;
    }


    public double probTrueVariant(double VQSRqual) {
        for ( VQSRRange r : points ) {
            if ( VQSRqual <= r.getStart() && VQSRqual > r.getStop() )
                return r.getTruePositiveRate();
        }

        throw new ReviewedGATKException("BUG: should not be able to reach this code");
    }

    public double probTrueVariant(String VQSRQualKey, VariantContext vc) {
        if ( vc.isFiltered() )
            return 0.0;
        else if ( vc.hasAttribute(VQSRQualKey) ) {
            double qual = vc.getAttributeAsDouble(VQSRQualKey, 0.0);
            return probTrueVariant(qual);
        } else {
            throw new UserException.VariantContextMissingRequiredField(VQSRQualKey, vc);
        }
    }

    /**
     * Returns a likelihoods vector adjusted by the probability that the site is an error.  Returns a
     * null vector if the probability of the site being real is 0.0
     * @param VQSRQualKey
     * @param vc
     * @param log10Likelihoods
     * @return
     */
    public double[] includeErrorRateInLikelihoods(String VQSRQualKey, VariantContext vc, double[] log10Likelihoods) {
        double[] updated = new double[log10Likelihoods.length];

        double alpha = probTrueVariant(VQSRQualKey, vc);

        if ( alpha == CERTAIN_FALSE_POSITIVE )
            return null;
        else {
            double noInfoPr = 1.0 / 3;
            if ( DEBUG ) System.out.printf("------------------------------%n");
            for ( int i = 0; i < log10Likelihoods.length; i++) {
                double p = Math.pow(10, log10Likelihoods[i]);
                double q = alpha * p + (1-alpha) * noInfoPr;
                if ( DEBUG ) System.out.printf("  vqslod = %.2f, p = %.2e, alpha = %.2e, q = %.2e%n", vc.getAttributeAsDouble(VQSRQualKey, 0.0), p, alpha, q);
                updated[i] = Math.log10(q);
            }

            return updated;
        }
    }


    public void printInfo(Logger logger) {
        for ( VQSRRange r : points ) {
            logger.info(String.format("  start=%f stop=%f TPrate=%.6e", r.getStart(), r.getStop(), r.getTruePositiveRate()));
        }
    }
}
