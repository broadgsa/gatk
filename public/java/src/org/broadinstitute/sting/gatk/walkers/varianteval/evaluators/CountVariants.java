package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.Analysis;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.DataPoint;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

@Analysis(description = "Counts different classes of variants in the sample")
public class CountVariants extends VariantEvaluator implements StandardEval {

    // the following fields are in output order:

    // basic counts on various rates found
    @DataPoint(description = "Number of processed loci")
    public long nProcessedLoci = 0;
    @DataPoint(description = "Number of called loci")
    public long nCalledLoci = 0;
    @DataPoint(description = "Number of reference loci")
    public long nRefLoci = 0;
    @DataPoint(description = "Number of variant loci")
    public long nVariantLoci = 0;

    // the following two calculations get set in the finalizeEvaluation
    @DataPoint(description = "Variants per loci rate")
    public double variantRate = 0;
    @DataPoint(description = "Number of variants per base")
    public double variantRatePerBp = 0;


    @DataPoint(description = "Number of snp loci")
    public long nSNPs = 0;
    @DataPoint(description = "Number of mnp loci")
    public long nMNPs = 0;
    @DataPoint(description = "Number of insertions")
    public long nInsertions = 0;
    @DataPoint(description = "Number of deletions")
    public long nDeletions = 0;
    @DataPoint(description = "Number of complex loci")
    public long nComplex = 0;


    @DataPoint(description = "Number of no calls loci")
    public long nNoCalls = 0;
    @DataPoint(description = "Number of het loci")
    public long nHets = 0;
    @DataPoint(description = "Number of hom ref loci")
    public long nHomRef = 0;
    @DataPoint(description = "Number of hom var loci")
    public long nHomVar = 0;
    @DataPoint(description = "Number of singletons")
    public long nSingletons = 0;
    @DataPoint(description = "Number of derived homozygotes")
    public long nHomDerived = 0;

    // calculations that get set in the finalizeEvaluation method
    @DataPoint(description = "heterozygosity per locus rate")
    public double heterozygosity = 0;
    @DataPoint(description = "heterozygosity per base pair")
    public double heterozygosityPerBp = 0;
    @DataPoint(description = "heterozygosity to homozygosity ratio")
    public double hetHomRatio = 0;
    @DataPoint(description = "indel rate (insertion count + deletion count)")
    public double indelRate = 0;
    @DataPoint(description = "indel rate per base pair")
    public double indelRatePerBp = 0;
    @DataPoint(description = "deletion to insertion ratio")
    public double deletionInsertionRatio = 0;
    
    private double perLocusRate(long n) {
        return rate(n, nProcessedLoci);
    }

    private long perLocusRInverseRate(long n) {
        return inverseRate(n, nProcessedLoci);
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
        nCalledLoci++;

        if (vc1.isVariant()) nVariantLoci++;
        switch (vc1.getType()) {
            case NO_VARIATION:
                nRefLoci++;
                break;
            case SNP:
                nSNPs++;
                if (vc1.getAttributeAsBoolean("ISSINGLETON")) nSingletons++;
                break;
            case MNP:
                nMNPs++;
                if (vc1.getAttributeAsBoolean("ISSINGLETON")) nSingletons++;
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

        String refStr = vc1.getReference().getBaseString().toUpperCase();

        String aaStr = vc1.hasAttribute("ANCESTRALALLELE") ? vc1.getAttributeAsString("ANCESTRALALLELE").toUpperCase() : null;
//        if (aaStr.equals(".")) {
//            aaStr = refStr;
//        }

        // ref  aa  alt  class
        // A    C   A    der homozygote
        // A    C   C    anc homozygote

        // A    A   A    ref homozygote
        // A    A   C
        // A    C   A
        // A    C   C

        for (Genotype g : vc1.getGenotypes().values()) {
            String altStr = vc1.getAlternateAlleles().size() > 0 ? vc1.getAlternateAllele(0).getBaseString().toUpperCase() : null;

            switch (g.getType()) {
                case NO_CALL:
                    nNoCalls++;
                    break;
                case HOM_REF:
                    nHomRef++;

                    if ( aaStr != null && altStr != null && !refStr.equalsIgnoreCase(aaStr) ) {
                        nHomDerived++;
                    }

                    break;
                case HET:
                    nHets++;
                    break;
                case HOM_VAR:
                    nHomVar++;

                    if ( aaStr != null && altStr != null && !altStr.equalsIgnoreCase(aaStr) ) {
                        nHomDerived++;
                    }

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