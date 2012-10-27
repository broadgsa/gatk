package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.GenotypeAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFormatHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFStandardHeaderLines;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypeBuilder;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * The depth of coverage of each VCF allele in this sample.
 *
 * The AD and DP are complementary fields that are two important ways of thinking about the depth of the data for this
 * sample at this site.  While the sample-level (FORMAT) DP field describes the total depth of reads that passed the
 * Unified Genotyper's internal quality control metrics (like MAPQ > 17, for example), the AD values (one for each of
 * REF and ALT fields) is the unfiltered count of all reads that carried with them the
 * REF and ALT alleles. The reason for this distinction is that the DP is in some sense reflective of the
 * power I have to determine the genotype of the sample at this site, while the AD tells me how many times
 * I saw each of the REF and ALT alleles in the reads, free of any bias potentially introduced by filtering
 * the reads. If, for example, I believe there really is a an A/T polymorphism at a site, then I would like
 * to know the counts of A and T bases in this sample, even for reads with poor mapping quality that would
 * normally be excluded from the statistical calculations going into GQ and QUAL. Please note, however, that
 * the AD isn't necessarily calculated exactly for indels. Only reads which are statistically favoring one allele over the other are counted.
 * Because of this fact, the sum of AD may be different than the individual sample depth, especially when there are
 * many non-informatice reads.
 * Because the AD includes reads and bases that were filtered by the Unified Genotyper and in case of indels is based on a statistical computation,
 * <b>one should not base assumptions about the underlying genotype based on it</b>;
 * instead, the genotype likelihoods (PLs) are what determine the genotype calls.
 */
public class DepthPerAlleleBySample extends GenotypeAnnotation implements StandardAnnotation {

    public void annotate(final RefMetaDataTracker tracker,
                         final AnnotatorCompatible walker,
                         final ReferenceContext ref,
                         final AlignmentContext stratifiedContext,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final PerReadAlleleLikelihoodMap alleleLikelihoodMap) {
        if ( g == null || !g.isCalled() || ( stratifiedContext == null && alleleLikelihoodMap == null) )
            return;

        if (alleleLikelihoodMap != null && !alleleLikelihoodMap.isEmpty())
            annotateWithLikelihoods(alleleLikelihoodMap, vc, gb);
        else if ( stratifiedContext != null && (vc.isSNP()))
            annotateWithPileup(stratifiedContext, vc, gb);
    }

    private void annotateWithPileup(final AlignmentContext stratifiedContext, final VariantContext vc, final GenotypeBuilder gb) {

        HashMap<Byte, Integer> alleleCounts = new HashMap<Byte, Integer>();
        for ( Allele allele : vc.getAlleles() )
            alleleCounts.put(allele.getBases()[0], 0);

        ReadBackedPileup pileup = stratifiedContext.getBasePileup();
        for ( PileupElement p : pileup ) {
            if ( alleleCounts.containsKey(p.getBase()) )
                alleleCounts.put(p.getBase(), alleleCounts.get(p.getBase())+p.getRepresentativeCount());
        }

        // we need to add counts in the correct order
        int[] counts = new int[alleleCounts.size()];
        counts[0] = alleleCounts.get(vc.getReference().getBases()[0]);
        for (int i = 0; i < vc.getAlternateAlleles().size(); i++)
            counts[i+1] = alleleCounts.get(vc.getAlternateAllele(i).getBases()[0]);

        gb.AD(counts);
    }

    private void annotateWithLikelihoods(final PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap, final VariantContext vc, final GenotypeBuilder gb) {
        final HashMap<Allele, Integer> alleleCounts = new HashMap<Allele, Integer>();

        for ( final Allele allele : vc.getAlleles() ) {
            alleleCounts.put(allele, 0);
        }
        for (Map.Entry<GATKSAMRecord,Map<Allele,Double>> el : perReadAlleleLikelihoodMap.getLikelihoodReadMap().entrySet()) {
            final GATKSAMRecord read = el.getKey();
            final Allele a = PerReadAlleleLikelihoodMap.getMostLikelyAllele(el.getValue());
            if (a.isNoCall())
                continue; // read is non-informative
            if (!vc.getAlleles().contains(a))
                continue; // sanity check - shouldn't be needed
            alleleCounts.put(a, alleleCounts.get(a) + (read.isReducedRead() ? read.getReducedCount(ReadUtils.getReadCoordinateForReferenceCoordinate(read, vc.getStart(), ReadUtils.ClippingTail.RIGHT_TAIL)) : 1));
        }
        final int[] counts = new int[alleleCounts.size()];
        counts[0] = alleleCounts.get(vc.getReference());
        for (int i = 0; i < vc.getAlternateAlleles().size(); i++)
            counts[i+1] = alleleCounts.get( vc.getAlternateAllele(i) );

        gb.AD(counts);
    }

    public List<String> getKeyNames() { return Arrays.asList(VCFConstants.GENOTYPE_ALLELE_DEPTHS); }

    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(VCFStandardHeaderLines.getFormatLine(getKeyNames().get(0)));
    }
}