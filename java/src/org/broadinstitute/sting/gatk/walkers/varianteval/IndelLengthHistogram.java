package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.playground.utils.report.tags.DataPoint;
import org.broadinstitute.sting.playground.utils.report.utils.TableType;
import org.broadinstitute.sting.utils.StingException;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date May 26, 2010
 */
public class IndelLengthHistogram extends VariantEvaluator {
    private final int SIZE_LIMIT = 50;
    @DataPoint(name="indelLengthHistogram",description="Histogram of indel lengths")
    IndelHistogram indelHistogram = new IndelHistogram(SIZE_LIMIT);

    /*
     * Indel length histogram table object
     */

    class IndelHistogram implements TableType {
        private Integer[] colKeys;
        private int limit;
        private String[] rowKeys = {"EventLength"};
        private int[] indelHistogram;

        public IndelHistogram(int limit) {
            colKeys = initColKeys(limit);
            indelHistogram = new int[limit*2];
            this.limit = limit;
        }

        public Object[] getColumnKeys() {
            return colKeys;
        }

        public Object[] getRowKeys() {
            return rowKeys;
        }

        public Object getCell(int row, int col) {
            return indelHistogram[col];
        }

        private Integer[] initColKeys(int size) {
            Integer[] cK = new Integer[size*2+1];
            int index = 0;
            for ( int i = -size; i < size; i ++ ) {
                cK[index] = i;
                index++;
            }

            return cK;
        }

        public String getName() { return "indelHistTable"; }

        public void update(int eLength) {
            indelHistogram[len2index(eLength)]++;
        }

        private int len2index(int len) {
            if ( len > limit || len < -limit ) {
                throw new StingException("Indel length exceeds limit of "+limit+" please increase indel limit size");
            }
            return len + limit;
        }
    }

    public IndelLengthHistogram(VariantEvalWalker parent) { super(parent); }

    public boolean enabled() { return false; }

    public String getName() { return "IndelLengthHistogram"; }

    public int getComparisonOrder() { return 1; } // need only the evals

    public String update1(VariantContext vc1, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( ! vc1.isBiallelic() && vc1.isIndel() ) {
            veWalker.getLogger().warn("[IndelLengthHistogram] Non-biallelic indel at "+ref.getLocus()+" ignored.");
            return vc1.toString(); // biallelic sites are output
        }

        if ( vc1.isIndel() ) {
            if ( vc1.isInsertion() ) {
                indelHistogram.update(vc1.getAlternateAllele(0).length());
            } else if ( vc1.isDeletion() ) {
                indelHistogram.update(-vc1.getReference().length());
            } else {
                throw new StingException("Indel type that is not insertion or deletion.");
            }
        }

        return null;
    }
}
