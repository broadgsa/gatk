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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.Analysis;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.DataPoint;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.TableType;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * @author rpoplin
 * @since Apr 6, 2010
 */

@Analysis(name = "Variant Quality Score", description = "Shows various stats of sets of variants binned by variant quality score")
public class VariantQualityScore extends VariantEvaluator {

    // a mapping from quality score histogram bin to Ti/Tv ratio
    @DataPoint(description = "the Ti/Tv ratio broken out by variant quality")
    TiTvStats titvStats = null;

    @DataPoint(description = "average variant quality for each allele count")
    AlleleCountStats alleleCountStats = null;

    static class TiTvStats implements TableType {
        final static int NUM_BINS = 20;
        final HashMap<Integer, Pair<Long,Long>> qualByIsTransition = new HashMap<Integer, Pair<Long,Long>>(); // A hashMap holds all the qualities until we are able to bin them appropriately
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
            StringBuffer returnString = new StringBuffer();
            // output the ti/tv array
            returnString.append("titvByQuality: ");
            for( int iii = 0; iii < NUM_BINS; iii++ ) {
                returnString.append(titvByQuality[iii]);
                returnString.append(" ");
            }
            return returnString.toString();
        }

        public void incrValue( final double qual, final boolean isTransition ) {
            final Integer qualKey = Math.round((float) qual);
            final long numTransition = (isTransition ? 1L : 0L);
            final long numTransversion = (isTransition ? 0L : 1L);
            if( qualByIsTransition.containsKey(qualKey) ) {
                Pair<Long,Long> transitionPair = qualByIsTransition.get(qualKey);
                transitionPair.set(transitionPair.getFirst() + numTransition, transitionPair.getSecond() + numTransversion);
                qualByIsTransition.put(qualKey, transitionPair);
            } else {
                qualByIsTransition.put(qualKey, new Pair<Long,Long>(numTransition,numTransversion));
            }
        }

        public void organizeTiTvTables() {
            for( int iii = 0; iii < NUM_BINS; iii++ ) {
                transitionByQuality[iii] = 0L;
                transversionByQuality[iii] = 0L;
                titvByQuality[iii] = 0.0;
            }

            int maxQual = 0;

            // Calculate the maximum quality score in order to normalize and histogram
            for( final Integer qual : qualByIsTransition.keySet() ) {
                if( qual > maxQual ) {
                    maxQual = qual;
                }
            }

            final double binSize = ((double)maxQual) / ((double) (NUM_BINS-1));

            for( final Integer qual : qualByIsTransition.keySet() ) {
                final int index = (int)Math.floor( ((double) qual) / binSize );
                if( index >= 0 ) { // BUGBUG: why is there overflow here?
                    Pair<Long,Long> transitionPair = qualByIsTransition.get(qual);
                    transitionByQuality[index] += transitionPair.getFirst();
                    transversionByQuality[index] += transitionPair.getSecond();
                }
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

    class AlleleCountStats implements TableType {
        final HashMap<Integer, ArrayList<Double>> qualityListMap = new HashMap<Integer, ArrayList<Double>>();
        final HashMap<Integer, Double> qualityMap = new HashMap<Integer, Double>();

        public Object[] getRowKeys() {
            final int NUM_BINS = qualityListMap.keySet().size();
            final String rowKeys[] = new String[NUM_BINS];
            int iii = 0;
            for( final Integer key : qualityListMap.keySet() ) {
                rowKeys[iii] = "AC" + key;
                iii++;
            }
            return rowKeys;

        }

        public Object[] getColumnKeys() {
            return new String[]{"alleleCount","avgQual"};
        }

        public String getName() {
            return "AlleleCountStats";
        }

        public String getCell(int x, int y) {
            int iii = 0;
            for( final Integer key : qualityListMap.keySet() ) {
                if(iii == x) {
                    if(y == 0) { return String.valueOf(key); }
                    else { return String.valueOf(qualityMap.get(key)); }
                }
                iii++;
            }
            return null;
        }

        public String toString() {
            String returnString = "";
            // output the quality map
            returnString += "AlleleCountStats: ";
            //for( int iii = 0; iii < NUM_BINS; iii++ ) {
            //    returnString += titvByQuality[iii] + " ";
            //}
            return returnString;
        }

        public void incrValue( final double qual, final int alleleCount ) {
            ArrayList<Double> list = qualityListMap.get(alleleCount);
            if(list==null) { list = new ArrayList<Double>(); }
            list.add(qual);
            qualityListMap.put(alleleCount, list);
        }

        public void organizeAlleleCountTables() {
            for( final Integer key : qualityListMap.keySet() ) {
                final ArrayList<Double> list = qualityListMap.get(key);
                double meanQual = 0.0;
                final double numQuals = (double)list.size();
                for( Double qual : list ) {
                    meanQual += qual / numQuals;
                }
                qualityMap.put(key, meanQual);
            }
        }
    }

    //public VariantQualityScore(VariantEvalWalker parent) {
        //super(parent);
    //}

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

        if( eval != null && eval.isSNP() && eval.isBiallelic() && eval.isPolymorphic() ) { //BUGBUG: only counting biallelic sites (revisit what to do with triallelic sites)
            if( titvStats == null ) { titvStats = new TiTvStats(); }
            titvStats.incrValue(eval.getPhredScaledQual(), VariantContextUtils.isTransition(eval));

            if( alleleCountStats == null ) { alleleCountStats = new AlleleCountStats(); }
            int alternateAlleleCount = 0;
            for (final Allele a : eval.getAlternateAlleles()) {
                alternateAlleleCount += eval.getChromosomeCount(a);
            }
            alleleCountStats.incrValue(eval.getPhredScaledQual(), alternateAlleleCount);
        }

        return interesting; // This module doesn't capture any interesting sites, so return null
    }

    public void finalizeEvaluation() {
        if( titvStats != null ) {
            titvStats.organizeTiTvTables();
        }
        if( alleleCountStats != null ) {
            alleleCountStats.organizeAlleleCountTables();
        }
    }
}