package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.report.tags.Analysis;
import org.broadinstitute.sting.utils.report.tags.DataPoint;
import org.broadinstitute.sting.utils.report.utils.TableType;

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

@Analysis(name = "Indel Statistics", description = "Shows various indel metrics and statistics")
public class IndelStatistics extends VariantEvaluator {
    @DataPoint(name="IndelStatistics", description = "Indel Statistics")
    IndelStats indelStats = null;


    
    private static final int INDEL_SIZE_LIMIT = 100;
    private static final int NUM_SCALAR_COLUMNS = 10;

    static int len2Index(int ind) {
        return ind+INDEL_SIZE_LIMIT+NUM_SCALAR_COLUMNS;
    }

    static int index2len(int ind) {
        return ind-INDEL_SIZE_LIMIT-NUM_SCALAR_COLUMNS;
    }
    static class IndelStats implements TableType {
        protected final static String ALL_SAMPLES_KEY = "allSamples";
        protected final static String[] COLUMN_KEYS;
        static {
            COLUMN_KEYS= new String[NUM_SCALAR_COLUMNS+2*INDEL_SIZE_LIMIT+1];
            COLUMN_KEYS[0] = "heterozygosity";
            COLUMN_KEYS[1] = "number_of_insertions";
            COLUMN_KEYS[2] = "number_of_deletions";
            COLUMN_KEYS[3] = "number_het_insertions";
            COLUMN_KEYS[4] = "number_homozygous_insertions";
            COLUMN_KEYS[5] = "number_het_deletions";
            COLUMN_KEYS[6] = "number_homozygous_deletions";
            COLUMN_KEYS[7] = "number of homozygous reference sites";
            COLUMN_KEYS[8] = "number of complex events";
            COLUMN_KEYS[9] = "number of long indels";

            for (int k=NUM_SCALAR_COLUMNS; k < NUM_SCALAR_COLUMNS+ 2*INDEL_SIZE_LIMIT+1; k++)
                COLUMN_KEYS[k] = "indel_size_len"+Integer.valueOf(index2len(k));
        }

        // map of sample to statistics
        protected final HashMap<String, double[]> indelSummary = new HashMap<String, double[]>();

        public IndelStats(final VariantContext vc) {
            indelSummary.put(ALL_SAMPLES_KEY, new double[COLUMN_KEYS.length]);
            for( final String sample : vc.getGenotypes().keySet() ) {
                indelSummary.put(sample, new double[COLUMN_KEYS.length]);
            }
        }

        /**
         *
         * @return one row per sample
         */
        public Object[] getRowKeys() {
            return indelSummary.keySet().toArray(new String[indelSummary.size()]);
        }
        public Object getCell(int x, int y) {
            final Object[] rowKeys = getRowKeys();
            return String.format("%4.2f",indelSummary.get(rowKeys[x])[y]);
        }

        /**
         * get the column keys
         * @return a list of objects, in this case strings, that are the column names
         */
        public Object[] getColumnKeys() {
            return COLUMN_KEYS;
        }

        public String getName() {
            return "IndelStats";
        }

        public int getComparisonOrder() {
            return 1;   // we only need to see each eval track
        }

         public String toString() {
            return getName();
        }
/*
        private double ratio(long numer, long denom) {
            return denom != 0L ? 100.0 * ( ((double)numer) / ((double)denom) ) : 0.0;
        }

        final long[] allSamplesNumerators = new long[COLUMN_KEYS.length];
        final long[] allSamplesDenominators = new long[COLUMN_KEYS.length];

        private void updateSummaries(int i, double[] summary, long numer, long denom ) {
            allSamplesNumerators[i] += numer;
            allSamplesDenominators[i] += denom;
            summary[i] = ratio(numer, denom);
        }
        */

        /*
         * increment the specified value
         */
        public void incrValue(VariantContext vc) {
            int eventLength = 0;
            boolean isInsertion = false, isDeletion = false;

            if ( vc.isInsertion() ) {
                eventLength = vc.getAlternateAllele(0).length();
                indelSummary.get(ALL_SAMPLES_KEY)[1]++;
                isInsertion = true;
            } else if ( vc.isDeletion() ) {
                indelSummary.get(ALL_SAMPLES_KEY)[2]++;
                eventLength = -vc.getReference().length();
                isDeletion = true;
            }
            else {
                indelSummary.get(ALL_SAMPLES_KEY)[8]++;
            }

            // make sure event doesn't overstep array boundaries
            if (Math.abs(eventLength) < INDEL_SIZE_LIMIT)
                indelSummary.get(ALL_SAMPLES_KEY)[len2Index(eventLength)]++;
            else
                indelSummary.get(ALL_SAMPLES_KEY)[9]++;


            for( final String sample : vc.getGenotypes().keySet() ) {
                if ( indelSummary.containsKey(sample) ) {
                    Genotype g = vc.getGenotype(sample);
                    boolean isVariant = (g.isCalled() && !g.isHomRef());
                    if (isVariant) {
                        // update ins/del count
                        if (isInsertion) {
                            indelSummary.get(sample)[1]++;
                        }
                        else if (isDeletion)
                            indelSummary.get(sample)[2]++;
                        else
                            indelSummary.get(sample)[8]++;

                        // update histogram
                        if (Math.abs(eventLength) < INDEL_SIZE_LIMIT)
                            indelSummary.get(sample)[len2Index(eventLength)]++;
                        else
                            indelSummary.get(sample)[9]++;

                        if (g.isHet())
                            if (isInsertion)
                                indelSummary.get(sample)[3]++;
                            else
                                indelSummary.get(sample)[5]++;
                        else
                            if (isInsertion)
                                indelSummary.get(sample)[4]++;
                            else
                                indelSummary.get(sample)[6]++;



                    }
                    else
                        indelSummary.get(sample)[7]++;
                }
            }


        }
    }

    public IndelStatistics(VariantEvalWalker parent) {
        super(parent);
        // don't do anything
    }

    public String getName() {
        return "IndelStatistics";
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
    public String update2(VariantContext eval, VariantContext validation, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, VariantEvalWalker.EvaluationContext group) {
        return null;
    }


    public String update1(VariantContext eval, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        if (eval != null ) {
            if ( indelStats == null ) {
                int nSamples = this.getVEWalker().getNSamplesForEval(eval);
                if ( nSamples != -1 )
                    indelStats = new IndelStats(eval);
            }

            if ( eval.isIndel() &&
                    eval.isBiallelic() &&
                    indelStats != null ) {
                indelStats.incrValue(eval);
            }
        }

        return null; // This module doesn't capture any interesting sites, so return null
    }
    public String update0(VariantContext eval, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        return null;
    }
    public void finalizeEvaluation() {
    //
       int k=0; 
    }

}
