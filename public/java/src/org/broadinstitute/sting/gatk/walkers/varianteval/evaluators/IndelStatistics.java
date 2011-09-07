package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.Analysis;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.DataPoint;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.TableType;
import org.broadinstitute.sting.utils.IndelUtils;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

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

@Analysis(name = "IndelStatistics", description = "Shows various indel metrics and statistics")
public class IndelStatistics extends VariantEvaluator {
    @DataPoint(description = "Indel Statistics")
    IndelStats indelStats = null;

   // @DataPoint(description = "Indel Classification")
    IndelClasses indelClasses = null;

    int numSamples = 0;

    public void initialize(VariantEvalWalker walker) {
        numSamples = walker.getNumSamples();
    }

    private static final int INDEL_SIZE_LIMIT = 100;
    private static final int IND_HET = 0;
    private static final int IND_INS = 1;
    private static final int IND_DEL = 2;
    private static final int IND_COMPLEX = 3;
    private static final int IND_HET_INS = 4;
    private static final int IND_HOM_INS = 5;
    private static final int IND_HET_DEL = 6;
    private static final int IND_HOM_DEL = 7;
    private static final int IND_HOM_REF = 8;
    private static final int IND_MIXED = 9;
    private static final int IND_LONG = 10;
    private static final int IND_AT_EXP = 11;
    private static final int IND_CG_EXP = 12;
    private static final int IND_FRAMESHIFT = 13;
    private static final int NUM_SCALAR_COLUMNS = 14;

    static int len2Index(int ind) {
        return ind+INDEL_SIZE_LIMIT+NUM_SCALAR_COLUMNS;
    }

    static int index2len(int ind) {
        return ind-INDEL_SIZE_LIMIT-NUM_SCALAR_COLUMNS;
    }

    static class IndelStats implements TableType {
         protected final static String[] COLUMN_KEYS;

         static {
            COLUMN_KEYS= new String[NUM_SCALAR_COLUMNS+2*INDEL_SIZE_LIMIT+1];
            COLUMN_KEYS[0] = "heterozygosity";
            COLUMN_KEYS[1] = "insertions";
            COLUMN_KEYS[2] = "deletions";
            COLUMN_KEYS[3] = "complex";
            COLUMN_KEYS[4] = "het_insertions";
            COLUMN_KEYS[5] = "homozygous_insertions";
            COLUMN_KEYS[6] = "het_deletions";
            COLUMN_KEYS[7] = "homozygous_deletions";
            COLUMN_KEYS[8] = "homozygous_reference_sites";
            COLUMN_KEYS[9] = "complex_events";
            COLUMN_KEYS[10] = "long_indels";
            COLUMN_KEYS[11] = "AT_expansions";
            COLUMN_KEYS[12] = "CG_expansions";
            COLUMN_KEYS[13] = "frameshift_indels";

            for (int k=NUM_SCALAR_COLUMNS; k < NUM_SCALAR_COLUMNS+ 2*INDEL_SIZE_LIMIT+1; k++)
                COLUMN_KEYS[k] = "indel_size_len"+Integer.valueOf(index2len(k));
        }

        // map of sample to statistics
        protected final  int[] indelSummary;

        public IndelStats(final VariantContext vc) {
            indelSummary = new int[COLUMN_KEYS.length];
        }

        /**
         *
         * @return one row per sample
         */
        public Object[] getRowKeys() {
            return new String[]{"all"};
        }
        public Object getCell(int x, int y) {
            return String.format("%d",indelSummary[y]);
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
        public void incrValue(VariantContext vc, ReferenceContext ref) {
            int eventLength = 0;
            boolean isInsertion = false, isDeletion = false;

            if ( vc.isSimpleInsertion() ) {
                eventLength = vc.getAlternateAllele(0).length();
                indelSummary[IND_INS]++;
                isInsertion = true;
            } else if ( vc.isSimpleDeletion() ) {
                indelSummary[IND_DEL]++;
                eventLength = -vc.getReference().length();
                isDeletion = true;
            }
            else if (vc.isComplexIndel()) {
                indelSummary[IND_COMPLEX]++;
            }
            else if (vc.isMixed())
                indelSummary[IND_MIXED]++;

            if (IndelUtils.isATExpansion(vc,ref))
                indelSummary[IND_AT_EXP]++;
            if (IndelUtils.isCGExpansion(vc,ref))
                 indelSummary[IND_CG_EXP]++;

            // make sure event doesn't overstep array boundaries
            if (vc.isSimpleDeletion() || vc.isSimpleInsertion()) {
                if (Math.abs(eventLength) < INDEL_SIZE_LIMIT) {
                    indelSummary[len2Index(eventLength)]++;
                    if (eventLength % 3 != 0)
                        indelSummary[IND_FRAMESHIFT]++;
                }
                else
                    indelSummary[IND_LONG]++;
            }

        }
    }

    static class IndelClasses implements TableType {
        protected final static String[] columnNames = IndelUtils.getIndelClassificationNames();


        // map of sample to statistics
        protected final int[] indelClassSummary;

        public IndelClasses(final VariantContext vc) {
            indelClassSummary = new int[columnNames.length];
        }

        /**
         *
         * @return one row per sample
         */
        public Object[] getRowKeys() {
            return new String[]{"all"};
        }
        public Object getCell(int x, int y) {
            return String.format("%d",indelClassSummary[y]);
        }

        /**
         * get the column keys
         * @return a list of objects, in this case strings, that are the column names
         */
        public Object[] getColumnKeys() {
            return columnNames;
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
            indelClassSummary[index]++;
        }
        /*
         * increment the specified value
         */
         public void incrValue(VariantContext vc, ReferenceContext ref) {


            ArrayList<Integer> indices = IndelUtils.findEventClassificationIndex(vc,ref);
             //System.out.format("pos:%d \nREF: %s, ALT: %s\n",vc.getStart(), vc.getReference().getDisplayString(),
             //  vc.getAlternateAllele(0).getDisplayString());

             byte[] refBases = ref.getBases();
             //System.out.format("ref bef:%s\n",new String(Arrays.copyOfRange(refBases,0,refBases.length/2+1) ));
             //System.out.format("ref aft:%s\n",new String(Arrays.copyOfRange(refBases,refBases.length/2+1,refBases.length) ));
            for (int index: indices)    {
                incrementSampleStat(vc, index);
               // System.out.println(IndelUtils.getIndelClassificationName(index));
            }
        }

    }

    //public IndelStatistics(VariantEvalWalker parent) {
        //super(parent);
        // don't do anything
    //}

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

    public String update1(VariantContext eval, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        if (eval != null && eval.isPolymorphic()) {
            if ( indelStats == null ) {
                indelStats = new IndelStats(eval);
            }
            if ( indelClasses == null ) {
                indelClasses = new IndelClasses(eval);
            }

            if ( eval.isIndel() || eval.isMixed() ) {
                if (indelStats != null )
                    indelStats.incrValue(eval, ref);

                if (indelClasses != null)
                    indelClasses.incrValue(eval, ref);
            }
        }

        return null; // This module doesn't capture any interesting sites, so return null
    }

    public void finalizeEvaluation() {
        int k=0;
    }

}
