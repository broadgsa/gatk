package org.broadinstitute.sting.oneoffprojects.walkers.varianteval2;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.playground.utils.report.tags.Analysis;
import org.broadinstitute.sting.playground.utils.report.tags.DataPoint;
import org.broadinstitute.sting.playground.utils.report.utils.TableType;

import java.util.ArrayList;
import java.util.HashMap;

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
 * @author rpoplin
 * @since Apr 6, 2010
 */

@Analysis(name = "Variant Quality Score", description = "Shows various stats of sets of variants binned by variant quality score")
public class VariantQualityScore extends VariantEvaluator {

    // a mapping from quality score histogram bin to Ti/Tv ratio
    @DataPoint(name="TiTv by Quality", description = "the Ti/Tv ratio broken out by variant quality")
    TiTvStats titvStats = null;

    //@DataPoint(name="Quality by Allele Count", description = "average variant quality for each allele count")
    //AlleleCountStats alleleCountStats = null;

    class TiTvStats implements TableType {
        final int NUM_BINS = 20;
        final ArrayList<Double> qualities = new ArrayList<Double>(); // An ArrayList holds all the qualities until we are able to bin them appropriately
        final ArrayList<Boolean> isTransition = new ArrayList<Boolean>();
        final long transitionByQuality[] = new long[NUM_BINS];
        final long transversionByQuality[] = new long[NUM_BINS];
        final double titvByQuality[] = new double[NUM_BINS]; // the final ti/tv sets that get reported out

        public Object[] getRowKeys() {
            return new String[]{"sample"};
        }

        public Object[] getColumnKeys() {
            final String columnKeys[] = new String[NUM_BINS];
            for( int iii = 0; iii < NUM_BINS; iii++ ) {
                columnKeys[iii] = "titvBin" + iii;
            }
            return columnKeys;
        }

        public String getName() {
            return "TiTvStats";
        }

        public String getCell(int x, int y) {
            return String.valueOf(titvByQuality[y]);
        }

        public String toString() {
            String returnString = "";
            // output the ti/tv array
            returnString += "titvByQuality: ";
            for( int iii = 0; iii < NUM_BINS; iii++ ) {
                returnString += titvByQuality[iii] + " ";
            }
            return returnString;
        }

        public void incrValue( final double qual, final boolean _isTransition ) {
            qualities.add(qual);
            isTransition.add(_isTransition);
        }

        public void organizeTiTvTables() {
            for( int iii = 0; iii < NUM_BINS; iii++ ) {
                transitionByQuality[iii] = 0L;
                transversionByQuality[iii] = 0L;
                titvByQuality[iii] = 0.0;
            }

            double maxQual = 0.0;

            // Calculate the maximum quality score in order to normalize and histogram
            for( final Double qual : qualities ) {
                if( qual > maxQual ) {
                    maxQual = qual;
                }
            }

            final double binSize = maxQual / ((double) (NUM_BINS-1));

            int jjj = 0;
            for( final Double qual : qualities ) {
                final int index = (int)Math.floor( qual / binSize );
                if(isTransition.get(jjj)) {
                    transitionByQuality[index]++;
                } else {
                    transversionByQuality[index]++;
                }
                jjj++;
            }

            for( int iii = 0; iii < NUM_BINS; iii++ ) {
                if( transitionByQuality[iii] + transversionByQuality[iii] > 800L ) { // need to have a sufficient number of variants to get a useful Ti/Tv ratio
                    titvByQuality[iii] = ((double) transitionByQuality[iii]) / ((double) transversionByQuality[iii]);
                } else {
                    titvByQuality[iii] = 0.0;
                }
            }

        }
    }

    /*
    class AlleleCountStats implements TableType {
        final HashMap<Integer, ArrayList<Double>> qualityListMap = new HashMap<Integer, ArrayList<Double>>();
        final HashMap<Integer, Double> qualityMap = new HashMap<Integer, Double>();

        public Object[] getRowKeys() {            
            return new String[]{"sample"};
        }

        public Object[] getColumnKeys() {
            final int NUM_BINS = qualityListMap.keySet().size();
            final String columnKeys[] = new String[NUM_BINS];
            int iii = 0;
            for( final Integer key : qualityListMap.keySet() ) {
                columnKeys[iii] = "AC" + key;
                iii++;
            }
            return columnKeys;
        }

        public String getName() {
            return "TiTvStats";
        }

        public String getCell(int x, int y) {
            return String.valueOf(titvByQuality[y]);
        }

        public String toString() {
            String returnString = "";
            // output the ti/tv array
            returnString += "titvByQuality: ";
            for( int iii = 0; iii < NUM_BINS; iii++ ) {
                returnString += titvByQuality[iii] + " ";
            }
            return returnString;
        }

        public void incrValue( final double qual, final boolean _isTransition ) {
            qualities.add(qual);
            isTransition.add(_isTransition);
        }

        public void organizeTiTvTables() {
            for( int iii = 0; iii < NUM_BINS; iii++ ) {
                transitionByQuality[iii] = 0L;
                transversionByQuality[iii] = 0L;
                titvByQuality[iii] = 0.0;
            }

            double maxQual = 0.0;

            // Calculate the maximum quality score in order to normalize and histogram
            for( final Double qual : qualities ) {
                if( qual > maxQual ) {
                    maxQual = qual;
                }
            }

            final double binSize = maxQual / ((double) (NUM_BINS-1));

            int jjj = 0;
            for( final Double qual : qualities ) {
                final int index = (int)Math.floor( qual / binSize );
                if(isTransition.get(jjj)) {
                    transitionByQuality[index]++;
                } else {
                    transversionByQuality[index]++;
                }
                jjj++;
            }

            for( int iii = 0; iii < NUM_BINS; iii++ ) {
                if( transitionByQuality[iii] + transversionByQuality[iii] > 800L ) { // need to have a sufficient number of variants to get a useful Ti/Tv ratio
                    titvByQuality[iii] = ((double) transitionByQuality[iii]) / ((double) transversionByQuality[iii]);
                } else {
                    titvByQuality[iii] = 0.0;
                }
            }

        }
    }
    */

    public VariantQualityScore(VariantEval2Walker parent) {
        // don't do anything
    }

    public String getName() {
        return "VariantQualityScore";
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

        if(eval != null) {
            if( titvStats == null ) { titvStats = new TiTvStats(); }
            titvStats.incrValue(eval.getPhredScaledQual(), eval.isTransition());
        }

        return interesting; // This module doesn't capture any interesting sites, so return null
    }

    public void finalizeEvaluation() {
        if( titvStats!=null ) {
            titvStats.organizeTiTvTables();
        }
    }
}