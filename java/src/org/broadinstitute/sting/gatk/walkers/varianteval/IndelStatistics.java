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

import java.util.Arrays;
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

@Analysis(name = "IndelStatistics", description = "Shows various indel metrics and statistics")
public class IndelStatistics extends VariantEvaluator {
    @DataPoint(name="IndelStatistics", description = "Indel Statistics")
    IndelStats indelStats = null;

    @DataPoint(name="IndelClasses", description = "Indel Classification")
    IndelClasses indelClasses = null;

    
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

    static class IndelClasses implements TableType {
        protected final static String ALL_SAMPLES_KEY = "allSamples";
        protected final static String[] COLUMN_KEYS;



        static {
            COLUMN_KEYS= new String[41];
            COLUMN_KEYS[0] = "Novel_A";
            COLUMN_KEYS[1] = "Novel_C";
            COLUMN_KEYS[2] = "Novel_G";
            COLUMN_KEYS[3] = "Novel_T";
            COLUMN_KEYS[4]  = "NOVEL_1";
            COLUMN_KEYS[5]  = "NOVEL_2";
            COLUMN_KEYS[6]  = "NOVEL_3";
            COLUMN_KEYS[7]  = "NOVEL_4";
            COLUMN_KEYS[8]  = "NOVEL_5";
            COLUMN_KEYS[9]  = "NOVEL_6";
            COLUMN_KEYS[10] = "NOVEL_7";
            COLUMN_KEYS[11] = "NOVEL_8";
            COLUMN_KEYS[12] = "NOVEL_9";
            COLUMN_KEYS[13] = "NOVEL_10orMore";
            COLUMN_KEYS[14] = "RepeatExpansion_A";
            COLUMN_KEYS[15] = "RepeatExpansion_C";
            COLUMN_KEYS[16] = "RepeatExpansion_G";
            COLUMN_KEYS[17] = "RepeatExpansion_T";
            COLUMN_KEYS[18] = "RepeatExpansion_AC";
            COLUMN_KEYS[19] = "RepeatExpansion_AG";
            COLUMN_KEYS[20] = "RepeatExpansion_AT";
            COLUMN_KEYS[21] = "RepeatExpansion_CA";
            COLUMN_KEYS[22] = "RepeatExpansion_CG";
            COLUMN_KEYS[23] = "RepeatExpansion_CT";
            COLUMN_KEYS[24] = "RepeatExpansion_GA";
            COLUMN_KEYS[25] = "RepeatExpansion_GC";
            COLUMN_KEYS[26] = "RepeatExpansion_GT";
            COLUMN_KEYS[27] = "RepeatExpansion_TA";
            COLUMN_KEYS[28] = "RepeatExpansion_TC";
            COLUMN_KEYS[29] = "RepeatExpansion_TG";
            COLUMN_KEYS[30] = "RepeatExpansion_1";
            COLUMN_KEYS[31] = "RepeatExpansion_2";
            COLUMN_KEYS[32] = "RepeatExpansion_3";
            COLUMN_KEYS[33] = "RepeatExpansion_4";
            COLUMN_KEYS[34] = "RepeatExpansion_5";
            COLUMN_KEYS[35] = "RepeatExpansion_6";
            COLUMN_KEYS[36] = "RepeatExpansion_7";
            COLUMN_KEYS[37] = "RepeatExpansion_8";
            COLUMN_KEYS[38] = "RepeatExpansion_9";
            COLUMN_KEYS[39] = "RepeatExpansion_10orMore";
            COLUMN_KEYS[40] = "Other";

        }

        private static final int START_IND_NOVEL = 4;
        private static final int STOP_IND_NOVEL = 13;
        private static final int START_IND_FOR_REPEAT_EXPANSION_1 = 14;
        private static final int STOP_IND_FOR_REPEAT_EXPANSION_2 = 29;
        private static final int START_IND_FOR_REPEAT_EXPANSION_COUNTS = 30;
        private static final int STOP_IND_FOR_REPEAT_EXPANSION_COUNTS = 39;
        private static final int IND_FOR_OTHER_EVENT = 40;
        private static final int START_IND_NOVEL_PER_BASE = 0;
        private static final int STOP_IND_NOVEL_PER_BASE = 3;


        // map of sample to statistics
        protected final HashMap<String, int[]> indelClassSummary = new HashMap<String, int[]>();

        public IndelClasses(final VariantContext vc) {
            indelClassSummary.put(ALL_SAMPLES_KEY, new int[COLUMN_KEYS.length]);
            for( final String sample : vc.getGenotypes().keySet() ) {
                indelClassSummary.put(sample, new int[COLUMN_KEYS.length]);
            }
        }

        /**
         *
         * @return one row per sample
         */
        public Object[] getRowKeys() {
            return indelClassSummary.keySet().toArray(new String[indelClassSummary.size()]);
        }
        public Object getCell(int x, int y) {
            final Object[] rowKeys = getRowKeys();
            return String.format("%d",indelClassSummary.get(rowKeys[x])[y]);
        }

        /**
         * get the column keys
         * @return a list of objects, in this case strings, that are the column names
         */
        public Object[] getColumnKeys() {
            return COLUMN_KEYS;
        }

        public String getName() {
            return "IndelClasses";
        }

        public int getComparisonOrder() {
            return 1;   // we only need to see each eval track
        }

         public String toString() {
            return getName();
        }

        private void incrementSampleStat(VariantContext vc, int index) {
            indelClassSummary.get(ALL_SAMPLES_KEY)[index]++;
            for( final String sample : vc.getGenotypes().keySet() ) {
                 if ( indelClassSummary.containsKey(sample) ) {
                     Genotype g = vc.getGenotype(sample);
                     boolean isVariant = (g.isCalled() && !g.isHomRef());
                     if (isVariant)
                         // update  count
                         indelClassSummary.get(sample)[index]++;

                 }
             }

        }
        /*
         * increment the specified value
         */
        public void incrValue(VariantContext vc, ReferenceContext ref) {
            int eventLength = 0;
            boolean isInsertion = false, isDeletion = false;
            String indelAlleleString;

            if ( vc.isInsertion() ) {
                eventLength = vc.getAlternateAllele(0).length();
                isInsertion = true;
                indelAlleleString = vc.getAlternateAllele(0).getDisplayString();
            } else if ( vc.isDeletion() ) {
                eventLength = vc.getReference().length();
                isDeletion = true;
                indelAlleleString = vc.getReference().getDisplayString();
            }
            else {
                incrementSampleStat(vc, IND_FOR_OTHER_EVENT);

                return;
            }

            byte[] refBases = ref.getBases();



            // See first if indel is a repetition of bases before current
            int indStart = refBases.length/2-eventLength+1;

            boolean done = false;
            int numRepetitions = 0;
            while (!done) {
                if (indStart < 0)
                    done = true;
                else {
                    String refPiece = new String(Arrays.copyOfRange(refBases,indStart,indStart+eventLength));
                    if (refPiece.matches(indelAlleleString))
                    {
                        numRepetitions++;
                        indStart = indStart - eventLength;
                    }
                    else
                        done = true;

                }
            }

            // now do it forward
            done = false;
            indStart = refBases.length/2+1;
            while (!done) {
                if (indStart + eventLength >= refBases.length)
                    break;
                else {
                    String refPiece = new String(Arrays.copyOfRange(refBases,indStart,indStart+eventLength));
                    if (refPiece.matches(indelAlleleString))
                    {
                        numRepetitions++;
                        indStart = indStart + eventLength;
                    }
                    else
                        done = true;

                }
            }

            if (numRepetitions == 0) {
                //unrepeated sequence from surroundings
                int ind = START_IND_NOVEL + (eventLength-1);
                if (ind > STOP_IND_NOVEL)
                    ind = STOP_IND_NOVEL;
                incrementSampleStat(vc, ind);

                if (eventLength == 1) {
                    // log single base indels additionally by base
                    String keyStr = "Novel_" + indelAlleleString;
                    int k;
                    for (k=START_IND_NOVEL_PER_BASE; k <= STOP_IND_NOVEL_PER_BASE; k++) {
                        if (keyStr.matches(COLUMN_KEYS[k]))
                            break;
                    }
                    // log now event
                    incrementSampleStat(vc, k);
                }
            }
            else {
                int ind = START_IND_FOR_REPEAT_EXPANSION_COUNTS + (numRepetitions-1);
                if (ind > STOP_IND_FOR_REPEAT_EXPANSION_COUNTS)
                    ind = STOP_IND_FOR_REPEAT_EXPANSION_COUNTS;
                    incrementSampleStat(vc, ind);

                if (eventLength<=2) {
                    // for single or dinucleotide indels, we further log the base in which they occurred
                    String keyStr = "RepeatExpansion_" + indelAlleleString;
                    int k;
                    for (k=START_IND_FOR_REPEAT_EXPANSION_1; k <= STOP_IND_FOR_REPEAT_EXPANSION_2; k++) {
                        if (keyStr.matches(COLUMN_KEYS[k]))
                            break;
                    }
                    // log now event
                    incrementSampleStat(vc, k);
                }

            }
 //g+
 /*
            System.out.format("RefBefore: %s\n",new String(refBefore));
            System.out.format("RefAfter: %s\n",new String(refAfter));
            System.out.format("Indel Allele: %s\n", indelAlleleString);
            System.out.format("Num Repetitions: %d\n", numRepetitions);
 */
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
            if ( indelClasses == null ) {
                indelClasses = new IndelClasses(eval);
            }

            if ( eval.isIndel() && eval.isBiallelic() ) {
                if (indelStats != null )
                    indelStats.incrValue(eval);

                if (indelClasses != null)
                    indelClasses.incrValue(eval, ref);
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
