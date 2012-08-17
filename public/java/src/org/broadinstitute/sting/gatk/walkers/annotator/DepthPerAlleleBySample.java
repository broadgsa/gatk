package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.GenotypeAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFormatHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFStandardHeaderLines;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypeBuilder;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;


/**
 * The depth of coverage of each VCF allele in this sample.
 *
 * This and DP are complementary fields that are two important ways of thinking about the depth of the data for this sample
 * at this site. The DP field describe the total depth of reads that passed the Unified Genotypers internal
 * quality control metrics (like MAPQ > 17, for example), whatever base was present in the read at this site.
 * The AD values (one for each of REF and ALT fields) is the count of all reads that carried with them the
 * REF and ALT alleles. The reason for this distinction is that the DP is in some sense reflective of the
 * power I have to determine the genotype of the sample at this site, while the AD tells me how many times
 * I saw each of the REF and ALT alleles in the reads, free of any bias potentially introduced by filtering
 * the reads. If, for example, I believe there really is a an A/T polymorphism at a site, then I would like
 * to know the counts of A and T bases in this sample, even for reads with poor mapping quality that would
 * normally be excluded from the statistical calculations going into GQ and QUAL. Please note, however, that
 * the AD isn't necessarily calculated exactly for indels (it counts as non-reference only those indels that
 * are actually present and correctly left-aligned in the alignments themselves). Because of this fact and
 * because the AD includes reads and bases that were filtered by the Unified Genotyper, <b>one should not base
 * assumptions about the underlying genotype based on it</b>; instead, the genotype likelihoods (PLs) are what
 * determine the genotype calls (see below).
 */
public class DepthPerAlleleBySample extends GenotypeAnnotation implements StandardAnnotation {

    public void annotate(final RefMetaDataTracker tracker,
                         final AnnotatorCompatible walker,
                         final ReferenceContext ref,
                         final AlignmentContext stratifiedContext,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb) {
        if ( g == null || !g.isCalled() )
            return;

        if ( vc.isSNP() )
            annotateSNP(stratifiedContext, vc, gb);
        else if ( vc.isIndel() )
            annotateIndel(stratifiedContext, ref.getBase(), vc, gb);
    }

    private void annotateSNP(final AlignmentContext stratifiedContext, final VariantContext vc, final GenotypeBuilder gb) {

        HashMap<Byte, Integer> alleleCounts = new HashMap<Byte, Integer>();
        for ( Allele allele : vc.getAlleles() )
            alleleCounts.put(allele.getBases()[0], 0);

        ReadBackedPileup pileup = stratifiedContext.getBasePileup();
        for ( PileupElement p : pileup ) {
            if ( alleleCounts.containsKey(p.getBase()) )
                alleleCounts.put(p.getBase(), alleleCounts.get(p.getBase())+1);
        }

        // we need to add counts in the correct order
        int[] counts = new int[alleleCounts.size()];
        counts[0] = alleleCounts.get(vc.getReference().getBases()[0]);
        for (int i = 0; i < vc.getAlternateAlleles().size(); i++)
            counts[i+1] = alleleCounts.get(vc.getAlternateAllele(i).getBases()[0]);

        gb.AD(counts);
    }

    private void annotateIndel(final AlignmentContext stratifiedContext, final byte refBase, final VariantContext vc, final GenotypeBuilder gb) {
        ReadBackedPileup pileup = stratifiedContext.getBasePileup();
        if ( pileup == null )
            return;

        final HashMap<Allele, Integer> alleleCounts = new HashMap<Allele, Integer>();
        final Allele refAllele = vc.getReference();

        for ( final Allele allele : vc.getAlleles() ) {
            alleleCounts.put(allele, 0);
        }

        for ( PileupElement p : pileup ) {
            if ( p.isBeforeInsertion() ) {

                final Allele insertion = Allele.create((char)refBase + p.getEventBases(), false);
                if ( alleleCounts.containsKey(insertion) ) {
                    alleleCounts.put(insertion, alleleCounts.get(insertion)+1);
                }

            } else if ( p.isBeforeDeletionStart() ) {
                if ( p.getEventLength() == refAllele.length() - 1 ) {
                    // this is indeed the deletion allele recorded in VC
                    final Allele deletion = Allele.create(refBase);
                    if ( alleleCounts.containsKey(deletion) ) {
                        alleleCounts.put(deletion, alleleCounts.get(deletion)+1);
                    }
                }
            } else if ( p.getRead().getAlignmentEnd() > vc.getStart() ) {
                alleleCounts.put(refAllele, alleleCounts.get(refAllele)+1);
            }
        }

        final int[] counts = new int[alleleCounts.size()];
        counts[0] = alleleCounts.get(refAllele);
        for (int i = 0; i < vc.getAlternateAlleles().size(); i++)
            counts[i+1] = alleleCounts.get( vc.getAlternateAllele(i) );

        gb.AD(counts);
    }

 //   public String getIndelBases()
    public List<String> getKeyNames() { return Arrays.asList(VCFConstants.GENOTYPE_ALLELE_DEPTHS); }

    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(VCFStandardHeaderLines.getFormatLine(getKeyNames().get(0)));
    }
}