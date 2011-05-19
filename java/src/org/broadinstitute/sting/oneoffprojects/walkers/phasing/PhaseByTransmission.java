package org.broadinstitute.sting.oneoffprojects.walkers.phasing;

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.util.*;

/**
 * Given genotypes for a trio, phases child by transmission.  Computes probability that the determined phase is correct
 * given that the genotypes for mom and dad are correct (useful if you want to use this to compare phasing accuracy, but
 * want to break that comparison down by phasing confidence in the truth set).  Optionally filters out sites where the
 * phasing is indeterminate.
 */
public class PhaseByTransmission extends RodWalker<Integer, Integer> {
    @Argument(shortName="f", fullName="familyPattern", required=true, doc="Pattern for the family structure (usage: mom+dad=child)")
    public String familyStr = null;

    @Argument(shortName="filter", fullName="filterPhaseIndeterminateSites", required=false, doc="Filters out sites that are phase-indeterminate")
    public Boolean filterPhaseIndeterminateSites = false;

    @Output
    protected VCFWriter vcfWriter = null;

    private String SAMPLE_NAME_MOM;
    private String SAMPLE_NAME_DAD;
    private String SAMPLE_NAME_CHILD;

    /**
     * Parse the familial relationship specification, and initialize VCF writer
     */
    public void initialize() {
        String[] pieces = familyStr.split("[\\+\\=]");

        SAMPLE_NAME_MOM = pieces[0];
        SAMPLE_NAME_DAD = pieces[1];
        SAMPLE_NAME_CHILD = pieces[2];

        Set<String> samples = new HashSet<String>();
        samples.add(SAMPLE_NAME_MOM);
        samples.add(SAMPLE_NAME_DAD);
        samples.add(SAMPLE_NAME_CHILD);

        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        headerLines.addAll(VCFUtils.getHeaderFields(this.getToolkit()));
        headerLines.add(new VCFInfoHeaderLine("TP", 1, VCFHeaderLineType.Float, "Probability that the phase is correct given that the genotypes are correct"));
        vcfWriter.writeHeader(new VCFHeader(headerLines, samples));
    }

    /**
     * Given genotypes for mom, dad, and child, determine the phase for the child by examining the possible transmitted alleles.
     *
     * @param mom    mom's genotype
     * @param dad    dad's genotype
     * @param child  child's genotype
     * @return  phased version of child's genotype, or null in the case that the phasing cannot be determined
     */
    private Genotype getPhasedChildGenotype(Genotype mom, Genotype dad, Genotype child) {
        List<Allele> momAlleles = mom.getAlleles();
        List<Allele> dadAlleles = dad.getAlleles();

        Set<Genotype> possiblePhasedChildGenotypes = new HashSet<Genotype>();

        double transmissionProb = QualityUtils.qualToProb(mom.getNegLog10PError())*QualityUtils.qualToProb(dad.getNegLog10PError());

        Map<String, Object> attributes = new HashMap<String, Object>();
        attributes.putAll(child.getAttributes());
        attributes.put("TP", transmissionProb);

        for (Allele momAllele : momAlleles) {
            for (Allele dadAllele : dadAlleles) {
                if (momAllele.isCalled() && dadAllele.isCalled()) {
                    List<Allele> possiblePhasedChildAlleles = new ArrayList<Allele>();
                    possiblePhasedChildAlleles.add(momAllele);
                    possiblePhasedChildAlleles.add(dadAllele);

                    Genotype possiblePhasedChildGenotype = new Genotype(child.getSampleName(), possiblePhasedChildAlleles, child.getNegLog10PError(), null, attributes, true);
                    possiblePhasedChildGenotypes.add(possiblePhasedChildGenotype);
                }
            }
        }

        for (Genotype possiblePhasedChildGenotype : possiblePhasedChildGenotypes) {
            if (child.sameGenotype(possiblePhasedChildGenotype, true)) {
                return possiblePhasedChildGenotype;
            }
        }

        return null;
    }

    /**
     * For each variant in the file, determine the phasing for the child and replace the child's genotype with the trio's genotype
     *
     * @param tracker  the reference meta-data tracker
     * @param ref      the reference context
     * @param context  the alignment context
     * @return null
     */
    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker != null) {
            Collection<VariantContext> vcs = tracker.getVariantContexts(ref, "variant", null, context.getLocation(), true, true);

            for (VariantContext vc : vcs) {
                Genotype mom = vc.getGenotype(SAMPLE_NAME_MOM);
                Genotype dad = vc.getGenotype(SAMPLE_NAME_DAD);
                Genotype child = vc.getGenotype(SAMPLE_NAME_CHILD);

                Genotype childPhased = getPhasedChildGenotype(mom, dad, child);

                Collection<Genotype> phasedGenotypes = new ArrayList<Genotype>();
                Collection<Allele> alleles = vc.getAlleles();
                Set<String> filters = new HashSet<String>();
                Map<String, Object> attributes = new HashMap<String, Object>();

                phasedGenotypes.add(mom);
                phasedGenotypes.add(dad);

                filters.addAll(vc.getFilters());
                attributes.putAll(vc.getAttributes());

                if (childPhased == null) {
                    phasedGenotypes.add(child);

                    if (filterPhaseIndeterminateSites) {
                        filters.add("PHASE_INDETERMINATE");
                    }
                } else {
                    phasedGenotypes.add(childPhased);
                }

                VariantContext newvc = new VariantContext("PhaseByTransmission", ref.getLocus().getContig(), vc.getStart(), vc.getStart(), alleles, phasedGenotypes, vc.getNegLog10PError(), filters, attributes);

                vcfWriter.add(newvc, ref.getBase());
            }
        }

        return null;
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    @Override
    public Integer reduceInit() {
        return null;
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    @Override
    public Integer reduce(Integer value, Integer sum) {
        return null;
    }
}
