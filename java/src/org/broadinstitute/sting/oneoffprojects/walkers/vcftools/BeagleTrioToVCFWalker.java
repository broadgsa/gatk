package org.broadinstitute.sting.oneoffprojects.walkers.vcftools;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.oneoffprojects.walkers.varianteval2.MendelianViolationEvaluator;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.vcf.VCFWriter;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * Test routine for new VariantContext object
 */
@Requires(value={DataSource.REFERENCE},referenceMetaData={@RMD(name="variants",type=ReferenceOrderedDatum.class), @RMD(name="beagle",type=BeagleROD.class)})
public class BeagleTrioToVCFWalker extends RodWalker<VariantContext, Long> {
    @Argument(shortName="trio", doc="If provide, treats the input VCF as a single record containing genotypes for a single trio; String formatted as dad+mom=child", required=false)
    protected String TRIO_STRUCTURE;

    @Argument(shortName="eth", fullName="excludeTripleHets", doc="If provide, sites that are triple hets calls will not be phased, regardless of Beagle's value", required=false)
    protected boolean dontPhaseTripleHets = false;

    int nTripletHets = 0;

    private MendelianViolationEvaluator.TrioStructure trio = null;

    private VCFWriter writer;
    private boolean headerWritten = false;
    private final static String TRACK_NAME = "variants";
    private final static String BEAGLE_NAME = "beagle";

    public void initialize() {
        trio = MendelianViolationEvaluator.parseTrioDescription(TRIO_STRUCTURE);
        writer = new VCFWriter(out);
    }

    public VariantContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        VariantContext vc = null;

        if ( ref != null ) {
            vc = tracker.getVariantContext(TRACK_NAME, null, context.getLocation(), false);
            BeagleROD beagle = tracker.lookup(BEAGLE_NAME,BeagleROD.class);

            if ( vc != null ) {
                if ( ! headerWritten ) {
                    RodVCF vcfrod = tracker.lookup(TRACK_NAME,RodVCF.class);
                    writer.writeHeader(vcfrod.getHeader());
                    headerWritten = true;
                }

                //System.out.printf("VCF: %s%n", tracker.lookup(TRACK_NAME, null));
                vc = maybePhaseVC(vc, beagle);
            }
        }

        if ( vc != null )
            writer.addRecord(VariantContextAdaptors.toVCF(vc, ref.getBase()));

        return vc;
    }

    private VariantContext maybePhaseVC(VariantContext unphased, BeagleROD beagle) {
        if ( beagle == null ) {
            return unphased;
        } else {
            Map<String, List<String>> bglData = beagle.getGenotypes();
            List<String> momBgl = bglData.get(trio.mom);
            List<String> dadBgl = bglData.get(trio.dad);

            Genotype unphasedMom = unphased.getGenotype(trio.mom);
            Genotype unphasedDad = unphased.getGenotype(trio.dad);
            Genotype unphasedKid = unphased.getGenotype(trio.child);

            if ( dontPhaseTripleHets && unphasedMom.isHet() && unphasedDad.isHet() && unphasedKid.isHet() ) {
		nTripletHets++;
                return unphased;
	    }
            else {
                Allele momTrans = unphased.getAllele(momBgl.get(0));
                Allele momUntrans = unphased.getAllele(momBgl.get(1));
                Allele dadTrans = unphased.getAllele(dadBgl.get(0));
                Allele dadUntrans = unphased.getAllele(dadBgl.get(1));

                Genotype momG = phaseGenotype(unphasedMom, Arrays.asList(momTrans, momUntrans));
                Genotype dadG = phaseGenotype(unphasedDad, Arrays.asList(dadTrans, dadUntrans));
                Genotype kidG = phaseGenotype(unphasedKid, Arrays.asList(momTrans, dadTrans));

                return new VariantContext(unphased.getName(), unphased.getLocation(), unphased.getAlleles(),
                        Arrays.asList(momG, dadG, kidG), unphased.getNegLog10PError(), unphased.getFilters(), unphased.getAttributes());
            }
        }
    }

    private Genotype phaseGenotype(Genotype base, List<Allele> alleles) {
        return new Genotype(base.getSampleName(), alleles, base.getNegLog10PError(), base.getFilters(), base.getAttributes(), true);
    }

    public Long reduceInit() {
        return 0L;
    }

    public Integer reduce(VariantContext point, Integer sum) {
        return sum;
    }

    public void onTraversalDone(Long result) {
        logger.info(String.format("Ignored phasing of %d het/het/het genotypes", nTripletHets));
        //logger.info(String.format("Converted %d (%.2f%%) of these sites", result.nConverted, (100.0 * result.nConverted) / result.nVariants));
    }

    public Long reduce(VariantContext vc, Long prevReduce) {
        return ( vc == null ? prevReduce : prevReduce+1);
    }
}
