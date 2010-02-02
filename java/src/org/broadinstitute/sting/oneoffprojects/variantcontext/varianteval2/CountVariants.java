package org.broadinstitute.sting.oneoffprojects.variantcontext.varianteval2;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.oneoffprojects.variantcontext.VariantContext;
import org.broadinstitute.sting.oneoffprojects.variantcontext.Genotype;
import org.broadinstitute.sting.utils.StingException;

import java.util.List;
import java.util.Arrays;

public class CountVariants extends VariantEvaluator {
    long nProcessedLoci = 0;
    long nCalledLoci = 0;
    long nVariantLoci = 0;
    long nRefLoci = 0;

    long nSNPs = 0;
    long nInsertions = 0;
    long nDeletions = 0;
    long nComplex = 0;

    long nNoCalls = 0;
    long nHets = 0;
    long nHomRef = 0;
    long nHomVar = 0;

    private double rate(long n) {
        return n / (1.0 * Math.max(nProcessedLoci, 1));
    }

    private long inverseRate(long n) {
        return n == 0 ? 0 : nProcessedLoci / Math.max(n, 1);
    }

    private double ratio(long num, long denom) {
        return ((double)num) / (Math.max(denom, 1));
    }

    public String getName() {
        return "counter";
    }

    public int getComparisonOrder() {
        return 1;   // we only need to see each eval track
    }

    public void update0(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) { 
        nProcessedLoci += context.getSkippedBases() + 1;
    }

    public void update1(VariantContext vc1, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        //nProcessedLoci++;
        nCalledLoci++;

        if ( vc1.isVariant() ) nVariantLoci++;
        switch ( vc1.getType() ) {
            case NO_VARIATION: nRefLoci++; break;
            case SNP: nSNPs++; break;
            case INDEL:
                if ( vc1.isInsertion() ) nInsertions++; else nDeletions++;
                break;
            case MIXED: nComplex++; break;
        }

        for ( Genotype g : vc1.getGenotypes().values()  ) {
            switch ( g.getType() ) {
                case NO_CALL: nNoCalls++; break;
                case HOM_REF: nHomRef++; break;
                case HET:     nHets++; break;
                case HOM_VAR: nHomVar++; break;
                default: throw new StingException("BUG: Unexpected genotype type: " + g);
            }
        }
    }

    public String toString() {
        return getName() + ": " + summaryLine();
    }

    private String summaryLine() {
        return String.format("%d %d %d %d " +
                "%.2e %d " +
                "%d %d %d %d " +
                "%d %d %d " +
                "%.2e %d %.2f " +
                "%.2f %d %.2f",
                nProcessedLoci, nCalledLoci, nRefLoci, nVariantLoci,
                rate(nVariantLoci), inverseRate(nVariantLoci),
                nSNPs, nDeletions, nInsertions, nComplex,
                nHomRef, nHets, nHomVar,
                rate(nHets), inverseRate(nHets), ratio(nHets, nHomVar),
                rate(nDeletions + nInsertions), inverseRate(nDeletions + nInsertions), ratio(nDeletions, nInsertions));
    }

    private static List<String> HEADER =
            Arrays.asList("nProcessedLoci", "nCalledLoci", "nRefLoci", "nVariantLoci",
                    "variantRate", "variantRatePerBp",
                    "nSNPs", "nDeletions", "nInsertions", "nComplex",
                    "nHomRefGenotypes", "nHetGenotypes", "nHomVarGenotypes",
                    "heterozygosity", "heterozygosityPerBp", "hetHomRatio",
                    "indelRate", "indelRatePerBp", "deletionInsertionRatio");

    // making it a table
    public List<String> getTableHeader() {
        return HEADER;
    }

    public List<List<String>> getTableRows() {
        return Arrays.asList(Arrays.asList(summaryLine().split(" ")));
    }
}