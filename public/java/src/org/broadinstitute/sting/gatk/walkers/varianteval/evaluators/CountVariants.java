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
    @DataPoint(description = "Number of complex indels")
    public long nComplex = 0;
    @DataPoint(description = "Number of mixed loci (loci that can't be classified as a SNP, Indel or MNP)")
    public long nMixed = 0;


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

        // Note from Eric:
        // This is really not correct.  What we really want here is a polymorphic vs. monomorphic count (i.e. on the Genotypes).
        // So in order to maintain consistency with the previous implementation (and the intention of the original author), I've
        // added in a proxy check for monomorphic status here.
        // Protect against case when vc only as no-calls too - can happen if we strafity by sample and sample as a single no-call.
       if ( vc1.isMonomorphic() ) {
            nRefLoci++;
        } else {
             switch (vc1.getType()) {
                case NO_VARIATION:
                    // shouldn't get here
                    break;
                case SNP:
                    nVariantLoci++;
                    nSNPs++;
                    if (vc1.getAttributeAsBoolean("ISSINGLETON", false)) nSingletons++;
                    break;
                case MNP:
                    nVariantLoci++;
                    nMNPs++;
                    if (vc1.getAttributeAsBoolean("ISSINGLETON", false)) nSingletons++;
                    break;
                case INDEL:
                    nVariantLoci++;
                    if (vc1.isSimpleInsertion())
                        nInsertions++;
                    else if (vc1.isSimpleDeletion())
                        nDeletions++;
                    else
                        nComplex++;
                    break;
                case MIXED:
                    nVariantLoci++;
                    nMixed++;
                    break;
                case SYMBOLIC:
                    // ignore symbolic alleles, but don't fail
                    // todo - consistent way of treating symbolic alleles thgoughout codebase?
                    break;
                default:
                    throw new ReviewedStingException("Unexpected VariantContext type " + vc1.getType());
            }
        }

        String refStr = vc1.getReference().getBaseString().toUpperCase();

        String aaStr = vc1.hasAttribute("ANCESTRALALLELE") ? vc1.getAttributeAsString("ANCESTRALALLELE", null).toUpperCase() : null;
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
        indelRate = perLocusRate(nDeletions + nInsertions + nComplex);
        indelRatePerBp = perLocusRInverseRate(nDeletions + nInsertions + nComplex);
        deletionInsertionRatio = ratio(nDeletions, nInsertions);
    }
}