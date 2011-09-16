package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.Analysis;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.DataPoint;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.TableType;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date May 26, 2010
 */
@Analysis(name = "Indel length histograms", description = "Shows the distrbution of insertion/deletion event lengths (negative for deletion, positive for insertion)")
public class IndelLengthHistogram extends VariantEvaluator {
    private static final int SIZE_LIMIT = 100;
    @DataPoint(description="Histogram of indel lengths")
    IndelHistogram indelHistogram = new IndelHistogram(SIZE_LIMIT);

    /*
     * Indel length histogram table object
     */

    static class IndelHistogram implements TableType {
        private Integer[] colKeys;
        private int limit;
        private String[] rowKeys = {"EventLength"};
        private Integer[] indelHistogram;

        public IndelHistogram(int limit) {
            colKeys = initColKeys(limit);
            indelHistogram = initHistogram(limit);
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
            for ( int i = -size; i <= size; i ++ ) {
                cK[index] = i;
                index++;
            }

            return cK;
        }

        private Integer[] initHistogram(int size) {
            Integer[] hist = new Integer[size*2+1];
            for ( int i = 0; i < 2*size+1; i ++ ) {
                hist[i] = 0;
            }

            return hist;
        }

        public String getName() { return "indelHistTable"; }

        public void update(int eLength) {
            indelHistogram[len2index(eLength)]++;
        }

        private int len2index(int len) {
            if ( len > limit || len < -limit ) {
                throw new ReviewedStingException("Indel length exceeds limit of "+limit+" please increase indel limit size");
            }
            return len + limit;
        }
    }

    public boolean enabled() { return true; }

    public String getName() { return "IndelLengthHistogram"; }

    public int getComparisonOrder() { return 1; } // need only the evals

    public String update1(VariantContext vc1, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        if ( vc1.isIndel() && vc1.isPolymorphic() ) {

            if ( ! vc1.isBiallelic() ) {
                //veWalker.getLogger().warn("[IndelLengthHistogram] Non-biallelic indel at "+ref.getLocus()+" ignored.");
                return vc1.toString(); // biallelic sites are output
            }

            // only count simple insertions/deletions, not complex indels
            if ( vc1.isSimpleInsertion() ) {
                indelHistogram.update(vc1.getAlternateAllele(0).length());
            } else if ( vc1.isSimpleDeletion() ) {
                indelHistogram.update(-vc1.getReference().length());
            }
        }

        return null;
    }
}
