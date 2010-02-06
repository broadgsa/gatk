package org.broadinstitute.sting.oneoffprojects.walkers.varianteval2;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
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

    public CountVariants(VariantEval2Walker parent) {
        // don't do anything
    }

    private double perLocusRate(long n)      { return rate(n, nProcessedLoci); }
    private long perLocusRInverseRate(long n) { return inverseRate(n, nProcessedLoci); }

    public String getName() {
        return "counter";
    }

    public boolean enabled() { return true; }

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
                perLocusRate(nVariantLoci), perLocusRInverseRate(nVariantLoci),
                nSNPs, nDeletions, nInsertions, nComplex,
                nHomRef, nHets, nHomVar,
                perLocusRate(nHets), perLocusRInverseRate(nHets), ratio(nHets, nHomVar),
                perLocusRate(nDeletions + nInsertions), perLocusRInverseRate(nDeletions + nInsertions), ratio(nDeletions, nInsertions));
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