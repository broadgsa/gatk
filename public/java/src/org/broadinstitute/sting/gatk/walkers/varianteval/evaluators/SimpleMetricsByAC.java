package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.Degeneracy;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.Sample;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.Analysis;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.DataPoint;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.StateKey;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.TableType;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

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
public class SimpleMetricsByAC extends VariantEvaluator implements StandardEval {
    // a mapping from quality score histogram bin to Ti/Tv ratio
    @DataPoint(description = "TiTv by allele count")
    MetricsByAc metrics = null;

    private final static Object[] METRIC_COLUMNS = {"AC", "nTi", "nTv", "n", "TiTv"};
    private int numSamples;

    class MetricsAtAC {
        public int ac = -1, nTi = 0, nTv = 0;

        public MetricsAtAC(int ac) { this.ac = ac; }

        public void update(VariantContext eval) {
            if ( VariantContextUtils.isTransition(eval) )
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

        public String getCell(int ac, int y) {
            return metrics.get(ac).getColumn(y);
        }

        public String toString() {
            return "";
        }

        public void incrValue( VariantContext eval ) {
            int ac = -1;
            
            if ( eval.hasGenotypes() )
                ac = eval.getChromosomeCount(eval.getAlternateAllele(0));
            else if ( eval.hasAttribute("AC") ) {
                ac = eval.getAttributeAsInt("AC", -1);
            }

            if ( ac != -1 ) {
                metrics.get(ac).update(eval);
            }
        }
    }

    public void initialize(VariantEvalWalker walker) {
        numSamples = walker.getNumSamples();
        metrics = new MetricsByAc(2*numSamples);
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
        if (numSamples == 0) {
            return null;
        }

        final String interesting = null;

        if (eval != null) {
            if ( metrics == null ) {
                int nSamples = numSamples;

                if ( nSamples != -1 ) {
                    metrics = new MetricsByAc(2 * nSamples);
                }
            }

            if ( eval.isSNP() && eval.isBiallelic() && eval.isPolymorphic() && metrics != null ) {
                metrics.incrValue(eval);
            }
        }

        return interesting; // This module doesn't capture any interesting sites, so return null
    }

    @Override
    public boolean stateIsApplicable(StateKey stateKey) {
        String sampleClassName = Sample.class.getSimpleName();
        String degeneracyClassName = Degeneracy.class.getSimpleName();

        //return !(stateKey.containsKey(sampleClassName) && !stateKey.get(sampleClassName).equalsIgnoreCase("all"));

        if (stateKey.containsKey(sampleClassName) && !stateKey.get(sampleClassName).equalsIgnoreCase("all")) {
            return false;
        }

        if (stateKey.containsKey(degeneracyClassName) && !stateKey.get(degeneracyClassName).equalsIgnoreCase("all")) {
            return false;
        }

        return true;
    }
}
