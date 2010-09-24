package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.report.tags.Analysis;
import org.broadinstitute.sting.utils.report.tags.DataPoint;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;


@Analysis(name = "Count Variants", description = "Counts different classes of variants in the sample")
public class CountVariants extends VariantEvaluator implements StandardEval {

    // the following fields are in output order:

    // basic counts on various rates found
    @DataPoint(description = "Number of processed loci")
    long nProcessedLoci = 0;
    @DataPoint(description = "Number of called loci")
    long nCalledLoci = 0;
    @DataPoint(description = "Number of reference loci")
    long nRefLoci = 0;
    @DataPoint(description = "Number of variant loci")
    long nVariantLoci = 0;

    // the following two calculations get set in the finalizeEvaluation
    @DataPoint(description = "Variants per loci rate")
    double variantRate = 0;
    @DataPoint(description = "Number of variants per base")
    double variantRatePerBp = 0;


    @DataPoint(description = "Number of snp loci")
    long nSNPs = 0;
    @DataPoint(description = "Number of insertions")
    long nInsertions = 0;
    @DataPoint(description = "Number of deletions")
    long nDeletions = 0;
    @DataPoint(description = "Number of complex loci")
    long nComplex = 0;

    @DataPoint(description = "Number of no calls loci")
    long nNoCalls = 0;
    @DataPoint(description = "Number of het loci")
    long nHets = 0;
    @DataPoint(description = "Number of hom ref loci")
    long nHomRef = 0;
    @DataPoint(description = "Number of hom var loci")
    long nHomVar = 0;

    // calculations that get set in the finalizeEvaluation method
    @DataPoint(description = "heterozygosity per locus rate")
    double heterozygosity = 0;
    @DataPoint(description = "heterozygosity per base pair")
    double heterozygosityPerBp = 0;
    @DataPoint(description = "heterozygosity to homozygosity ratio")
    double hetHomRatio = 0;
    @DataPoint(description = "indel rate (insertion count + deletion count)")
    double indelRate = 0;
    @DataPoint(description = "indel rate per base pair")
    double indelRatePerBp = 0;
    @DataPoint(description = "deletion to insertion ratio")
    double deletionInsertionRatio = 0;
    
    public CountVariants(VariantEvalWalker parent) {
        super(parent);
    }

    private double perLocusRate(long n) {
        return rate(n, nProcessedLoci);
    }

    private long perLocusRInverseRate(long n) {
        return inverseRate(n, nProcessedLoci);
    }

    public String getName() {
        return "counter";
    }

    public boolean enabled() {
        return true;
    }

    public int getComparisonOrder() {
        return 1;   // we only need to see each eval track
    }

    public void update0(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        nProcessedLoci += context.getSkippedBases() + (ref == null ? 0 : 1);
    }

    public String update1(VariantContext vc1, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        //nProcessedLoci++;
        nCalledLoci++;

        if (vc1.isVariant()) nVariantLoci++;
        switch (vc1.getType()) {
            case NO_VARIATION:
                nRefLoci++;
                break;
            case SNP:
                nSNPs++;
                break;
            case INDEL:
                if (vc1.isInsertion()) nInsertions++;
                else nDeletions++;
                break;
            case MIXED:
                nComplex++;
                break;
            default:
                throw new ReviewedStingException("Unexpected VariantContext type " + vc1.getType());
        }

        for (Genotype g : vc1.getGenotypes().values()) {
            switch (g.getType()) {
                case NO_CALL:
                    nNoCalls++;
                    break;
                case HOM_REF:
                    nHomRef++;
                    break;
                case HET:
                    nHets++;
                    break;
                case HOM_VAR:
                    nHomVar++;
                    break;
                default:
                    throw new ReviewedStingException("BUG: Unexpected genotype type: " + g);
            }
        }

        return null; // we don't capture any interesting sites
    }

    public void finalizeEvaluation() {
        variantRate = perLocusRate(nVariantLoci);
        variantRatePerBp = perLocusRInverseRate(nVariantLoci);
        heterozygosity = perLocusRate(nHets);
        heterozygosityPerBp = perLocusRInverseRate(nHets);
        hetHomRatio = ratio(nHets, nHomVar);
        indelRate = perLocusRate(nDeletions + nInsertions);
        indelRatePerBp = perLocusRInverseRate(nDeletions + nInsertions);
        deletionInsertionRatio = ratio(nDeletions, nInsertions);
    }
}