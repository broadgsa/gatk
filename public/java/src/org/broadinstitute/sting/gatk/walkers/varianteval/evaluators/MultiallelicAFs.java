/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.Analysis;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.DataPoint;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.TableType;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

@Analysis(description = "Evaluation summary for multi-allelic variants")
public class MultiallelicAFs extends VariantEvaluator {
    final protected static Logger logger = Logger.getLogger(MultiallelicAFs.class);

    public enum Type {
        SNP, INDEL
    }

    @DataPoint(description="Histogram of allele frequencies for most common SNP alternate allele")
    AFHistogram AFhistogramMaxSnp = new AFHistogram();

    @DataPoint(description="Histogram of allele frequencies for less common SNP alternate alleles")
    AFHistogram AFhistogramMinSnp = new AFHistogram();

    @DataPoint(description="Histogram of allele frequencies for most common Indel alternate allele")
    AFHistogram AFhistogramMaxIndel = new AFHistogram();

    @DataPoint(description="Histogram of allele frequencies for less common Indel alternate alleles")
    AFHistogram AFhistogramMinIndel = new AFHistogram();

    /*
     * AF histogram table object
     */
    static class AFHistogram implements TableType {
        private Object[] rowKeys, colKeys = {"count"};
        private int[] AFhistogram;

        private static final double AFincrement = 0.01;
        private static final int numBins = (int)(1.00 / AFincrement);

        public AFHistogram() {
            rowKeys = initRowKeys();
            AFhistogram = new int[rowKeys.length];
        }

        public Object[] getColumnKeys() {
            return colKeys;
        }

        public Object[] getRowKeys() {
            return rowKeys;
        }

        public Object getCell(int row, int col) {
            return AFhistogram[row];
        }

        private static Object[] initRowKeys() {
            ArrayList<String> keyList = new ArrayList<String>(numBins + 1);
            for ( double a = 0.00; a <= 1.01; a += AFincrement ) {
                keyList.add(String.format("%.2f", a));
            }
            return keyList.toArray();
        }

        public String getName() { return "AFHistTable"; }

        public void update(final double AF) {
            final int bin = (int)(numBins * MathUtils.round(AF, 2));
            AFhistogram[bin]++;
       }
    }

    public void initialize(VariantEvalWalker walker) {}

    @Override public boolean enabled() { return true; }

    public int getComparisonOrder() {
        return 2;
    }

    public void update0(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {}

    public String update2(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( eval == null || eval.isMonomorphicInSamples() )
            return null;

        if ( !eval.isBiallelic() )
            return null;

        // update counts
        switch ( eval.getType() ) {
            case SNP:
                updateAFhistogram(eval, AFhistogramMaxSnp, AFhistogramMinSnp);
                break;
            case INDEL:
                updateAFhistogram(eval, AFhistogramMaxIndel, AFhistogramMinIndel);
                break;
            default:
                throw new UserException.BadInput("Unexpected variant context type: " + eval);
        }

        return null; // we don't capture any interesting sites
    }

    private void updateAFhistogram(VariantContext vc, AFHistogram max, AFHistogram min) {

        final Object obj = vc.getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY, null);
        if ( obj == null || !(obj instanceof List) )
            return;

        List<String> list = (List<String>)obj;
        ArrayList<Double> AFs = new ArrayList<Double>(list.size());
        for ( String str : list ) {
            AFs.add(Double.valueOf(str));
        }

        Collections.sort(AFs);
        max.update(AFs.get(AFs.size()-1));
        for ( int i = 0; i < AFs.size() - 1; i++ )
            min.update(AFs.get(i));
    }
}
