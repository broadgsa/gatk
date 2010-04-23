package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.playground.utils.report.tags.Analysis;
import org.broadinstitute.sting.playground.utils.report.tags.DataPoint;
import org.broadinstitute.sting.playground.utils.report.utils.TableType;
import org.broadinstitute.sting.utils.StingException;
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
 * @author depristo
 * @since Apr 11, 2010
 */

@Analysis(name = "Quality Metrics by allele count", description = "Shows various stats binned by allele count")
public class SimpleMetricsByAC extends VariantEvaluator {
    // a mapping from quality score histogram bin to Ti/Tv ratio
    @DataPoint(name="TiTv by AC", description = "TiTv by allele count")
    MetricsByAc metrics = null;

    //@DataPoint(name="Quality by Allele Count", description = "average variant quality for each allele count")
    //AlleleCountStats alleleCountStats = null;

    private final static Object[] METRIC_COLUMNS = {"AC", "nTi", "nTv", "n", "Ti/Tv"};

    class MetricsAtAC {
        public int ac = -1, nTi = 0, nTv = 0;

        public MetricsAtAC(int ac) { this.ac = ac; }

        public void update(VariantContext eval) {
            if ( eval.isTransition() )
                nTi++;
            else
                nTv++;
        }

        // corresponding to METRIC_COLUMNS
        public String getColumn(int i) {
            switch (i) {
                case 0: return String.valueOf(ac);
                case 1: return String.valueOf(nTi);
                case 2: return String.valueOf(nTv);
                case 3: return String.valueOf(nTi + nTv);
                case 4: return String.valueOf(ratio(nTi, nTv));
                default:
                    throw new StingException("Unexpected column " + i);
            }
        }
    }

    class MetricsByAc implements TableType {
        int nchromosomes = -1;
        ArrayList<MetricsAtAC> metrics = new ArrayList<MetricsAtAC>();
        Object[] rows = null;

        public MetricsByAc( int nchromosomes ) {
            this.nchromosomes = nchromosomes;
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
            int an = eval.getChromosomeCount();
            if ( an == nchromosomes ) { // ignore sites with no calls
                int ac = eval.getChromosomeCount(eval.getAlternateAllele(0));
                metrics.get(ac).update(eval);
            }
        }
    }

    public SimpleMetricsByAC(VariantEvalWalker parent) {
        // don't do anything
    }

    public String getName() {
        return "SimpleMetricsByAC";
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

        if (eval != null && eval.isSNP() && eval.hasGenotypes() ) {
            if ( metrics == null )
                metrics = new MetricsByAc(2 * eval.getNSamples());
            metrics.incrValue(eval);
        }

        return interesting; // This module doesn't capture any interesting sites, so return null
    }

    //public void finalizeEvaluation() {
    //
    //}
}
