package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.GenotypeAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFormatHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineCount;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.pileup.ExtendedEventPileupElement;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedExtendedEventPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


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

    private static String REF_ALLELE = "REF";

    private static String DEL = "DEL"; // constant, for speed: no need to create a key string for deletion allele every time

    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker, ReferenceContext ref, AlignmentContext stratifiedContext, VariantContext vc, Genotype g) {
        if ( g == null || !g.isCalled() )
            return null;

        if ( vc.isSNP() )
            return annotateSNP(stratifiedContext, vc);
        if ( vc.isIndel() )
            return annotateIndel(stratifiedContext, vc);

        return null;
    }

    private Map<String,Object> annotateSNP(AlignmentContext stratifiedContext, VariantContext vc) {

        if ( ! stratifiedContext.hasBasePileup() ) return null;

        HashMap<Byte, Integer> alleleCounts = new HashMap<Byte, Integer>();
        for ( Allele allele : vc.getAlleles() )
            alleleCounts.put(allele.getBases()[0], 0);

        ReadBackedPileup pileup = stratifiedContext.getBasePileup();
        for ( PileupElement p : pileup ) {
            if ( alleleCounts.containsKey(p.getBase()) )
                alleleCounts.put(p.getBase(), alleleCounts.get(p.getBase())+1);
        }

        // we need to add counts in the correct order
        Integer[] counts = new Integer[alleleCounts.size()];
        counts[0] = alleleCounts.get(vc.getReference().getBases()[0]);
        for (int i = 0; i < vc.getAlternateAlleles().size(); i++)
            counts[i+1] = alleleCounts.get(vc.getAlternateAllele(i).getBases()[0]);

        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), counts);
        return map;
    }

    private Map<String,Object> annotateIndel(AlignmentContext stratifiedContext, VariantContext vc) {

        if ( ! stratifiedContext.hasExtendedEventPileup() ) {
            return null;
        }

        ReadBackedExtendedEventPileup pileup = stratifiedContext.getExtendedEventPileup();
        if ( pileup == null )
            return null;

        HashMap<String, Integer> alleleCounts = new HashMap<String, Integer>();
        alleleCounts.put(REF_ALLELE,0);
        Allele refAllele = vc.getReference();

        for ( Allele allele : vc.getAlternateAlleles() ) {

            if ( allele.isNoCall() ) {
                continue; // this does not look so good, should we die???
            }

            alleleCounts.put(getAlleleRepresentation(allele), 0);
        }

        for ( ExtendedEventPileupElement e : pileup.toExtendedIterable() ) {
            if ( e.isInsertion() ) {

                final String b =  e.getEventBases();
                if ( alleleCounts.containsKey(b) ) {
                    alleleCounts.put(b, alleleCounts.get(b)+1);
                }

            } else {
                if ( e.isDeletion() ) {
                    if ( e.getEventLength() == refAllele.length() ) {
                        // this is indeed the deletion allele recorded in VC
                        final String b = DEL;
                        if ( alleleCounts.containsKey(b) ) {
                            alleleCounts.put(b, alleleCounts.get(b)+1);
                        }
                    }
//                    else {
//                        System.out.print("   deletion of WRONG length found");
//                    }
                }
                else {
                    if ( e.getRead().getAlignmentEnd() <= vc.getStart() ) {
                        continue;
                    }
                    alleleCounts.put(REF_ALLELE,alleleCounts.get(REF_ALLELE)+1);
                }
            }
        }

        Integer[] counts = new Integer[alleleCounts.size()];
        counts[0] = alleleCounts.get(REF_ALLELE);
        for (int i = 0; i < vc.getAlternateAlleles().size(); i++)
            counts[i+1] = alleleCounts.get( getAlleleRepresentation(vc.getAlternateAllele(i)) );

        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), counts);

        //map.put(getKeyNames().get(0), counts);
        return map;
    }

    private String getAlleleRepresentation(Allele allele) {
        if ( allele.isNull() ) { // deletion wrt the ref
             return DEL;
        } else { // insertion, pass actual bases
            return allele.getBaseString();
        }

    }

 //   public String getIndelBases()
    public List<String> getKeyNames() { return Arrays.asList("AD"); }

    public List<VCFFormatHeaderLine> getDescriptions() { return Arrays.asList(new VCFFormatHeaderLine(getKeyNames().get(0), VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Allelic depths for the ref and alt alleles in the order listed")); }
}