package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.report.tags.Analysis;
import org.broadinstitute.sting.utils.report.tags.DataPoint;
import org.broadinstitute.sting.utils.report.utils.TableType;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.ArrayList;

/*
 * Copyright (c) 2010 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @author delangel
 * @since Apr 11, 2010
 */

@Analysis(name = "Indel Metrics by allele count", description = "Shows various stats binned by allele count")
public class IndelMetricsByAC extends VariantEvaluator implements StandardEval {
    // a mapping from quality score histogram bin to Ti/Tv ratio
    @DataPoint(name="Indel Metrics by AC", description = "Indel Metrics by allele count")
    MetricsByAc metrics = null;

    //@DataPoint(name="Quality by Allele Count", description = "average variant quality for each allele count")
    //AlleleCountStats alleleCountStats = null;
    private static final int INDEL_SIZE_LIMIT = 100;
    private static final int NUM_SCALAR_COLUMNS = 6;
    static int len2Index(int ind) {
        return ind+INDEL_SIZE_LIMIT+NUM_SCALAR_COLUMNS;
    }

    static int index2len(int ind) {
        return ind-INDEL_SIZE_LIMIT-NUM_SCALAR_COLUMNS;
    }

    protected final static String[] METRIC_COLUMNS;
    static {
        METRIC_COLUMNS= new String[NUM_SCALAR_COLUMNS+2*INDEL_SIZE_LIMIT+1];
        METRIC_COLUMNS[0] = "AC";
        METRIC_COLUMNS[1] = "nIns";
        METRIC_COLUMNS[2] = "nDels";
        METRIC_COLUMNS[3] = "n";
        METRIC_COLUMNS[4] = "nComplex";
        METRIC_COLUMNS[5] = "nLong";

        for (int k=NUM_SCALAR_COLUMNS; k < NUM_SCALAR_COLUMNS+ 2*INDEL_SIZE_LIMIT+1; k++)
            METRIC_COLUMNS[k] = "indel_size_len"+Integer.valueOf(index2len(k));
    }

    class MetricsAtAC {
        public int ac = -1, nIns =0, nDel = 0, nComplex = 0, nLong;
        public int sizeCount[] = new int[2*INDEL_SIZE_LIMIT+1];

        public MetricsAtAC(int ac) { this.ac = ac; }

        public void update(VariantContext eval) {
            int eventLength = 0;
            if ( eval.isInsertion() ) {
                eventLength = eval.getAlternateAllele(0).length();
                nIns++;
            } else if ( eval.isDeletion() ) {
                eventLength = -eval.getReference().length();
                nDel++;
            }
            else {
                nComplex++;
            }
            if (Math.abs(eventLength) < INDEL_SIZE_LIMIT)
                sizeCount[len2Index(eventLength)]++;
            else
                nLong++;



         }

        // corresponding to METRIC_COLUMNS
        public String getColumn(int i) {
            if (i >= NUM_SCALAR_COLUMNS && i <=NUM_SCALAR_COLUMNS+ 2*INDEL_SIZE_LIMIT)
                return String.valueOf(sizeCount[i-NUM_SCALAR_COLUMNS]);

            switch (i) {
                case 0: return String.valueOf(ac);
                case 1: return String.valueOf(nIns);
                case 2: return String.valueOf(nDel);
                case 3: return String.valueOf(nIns + nDel);
                case 4: return String.valueOf(nComplex);
                case 5: return String.valueOf(nLong);

                default:
                    throw new ReviewedStingException("Unexpected column " + i);
            }
        }
    }

    class MetricsByAc implements TableType {
        ArrayList<MetricsAtAC> metrics = new ArrayList<MetricsAtAC>();
        Object[] rows = null;

        public MetricsByAc( int nchromosomes ) {
            rows = new Object[nchromosomes+1];
            metrics = new ArrayList<MetricsAtAC>(nchromosomes+1);
            for ( int i = 0; i < nchromosomes + 1; i++ ) {
                metrics.add(new MetricsAtAC(i));
                rows[i] = "ac" + i;
            }
        }

        public Object[] getRowKeys() {
            return rows;
        }

        public Object[] getColumnKeys() {
            return METRIC_COLUMNS;
        }

        public String getName() {
            return "MetricsByAc";
        }

        //
        public String getCell(int ac, int y) {
            return metrics.get(ac).getColumn(y);
        }

        public String toString() {
            String returnString = "";
            return returnString;
        }

        public void incrValue( VariantContext eval ) {
            int ac = -1;

            if ( eval.hasGenotypes() )
                ac = eval.getChromosomeCount(eval.getAlternateAllele(0));
            else if ( eval.hasAttribute("AC") ) {
                ac = Integer.valueOf(eval.getAttributeAsString("AC"));
            }

            if ( ac != -1 )
                metrics.get(ac).update(eval);
        }
    }

    public IndelMetricsByAC(VariantEvalWalker parent) {
        super(parent);
        // don't do anything
    }

    public String getName() {
        return "IndelMetricsByAC";
    }

    public int getComparisonOrder() {
        return 1;   // we only need to see each eval track
    }

    public boolean enabled() {
        return true;
    }

    public String toString() {
        return getName();
    }

    public String update1(VariantContext eval, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        final String interesting = null;

        if (eval != null ) {
            if ( metrics == null ) {
                int nSamples = this.getVEWalker().getNSamplesForEval(eval);
                if ( nSamples != -1 )
                    metrics = new MetricsByAc(2 * nSamples);
            }

            if ( eval.isIndel() && eval.isBiallelic() &&
                    metrics != null ) {
                metrics.incrValue(eval);
            }
        }

        return interesting; // This module doesn't capture any interesting sites, so return null
    }

    //public void finalizeEvaluation() {
    //
    //}
}
