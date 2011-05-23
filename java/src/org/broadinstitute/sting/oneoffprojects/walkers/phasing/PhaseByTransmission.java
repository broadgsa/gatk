package org.broadinstitute.sting.oneoffprojects.walkers.phasing;

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.util.*;

/**
 * Given genotypes for a trio, phases child by transmission.  Computes probability that the determined phase is correct
 * given that the genotypes for mom and dad are correct (useful if you want to use this to compare phasing accuracy, but
 * want to break that comparison down by phasing confidence in the truth set).  Optionally filters out sites where the
 * phasing is indeterminate.  This walker assumes there are only three samples in the VCF file to begin with.
 */
public class PhaseByTransmission extends RodWalker<Integer, Integer> {
    @Argument(shortName="f", fullName="familyPattern", required=true, doc="Pattern for the family structure (usage: mom+dad=child)")
    public String familyStr = null;

    @Argument(shortName="nofilters", fullName="disableFilters", required=false, doc="Disable filters for sites where the phase can't be determined, where the parental origin of the alleles is ambiguous (i.e. everyone is heterozygous), or Mendelian violations")
    public Boolean noFilters = false;

    @Output
    protected VCFWriter vcfWriter = null;

    private String SAMPLE_NAME_MOM;
    private String SAMPLE_NAME_DAD;
    private String SAMPLE_NAME_CHILD;

    private class AmbiguousAlleleOriginException extends Exception {}
    private class InsufficientInfoToPhaseGenotypeException extends Exception {}
    private class MendelianViolationException extends Exception {}

    private final String ROD_NAME = "variant";
    private final String AMBIGUOUS_ALLELE_ORIGIN_FILTER_NAME = "AmbiguousAlleleOrigin";
    private final String PHASE_INDETERMINATE_FILTER_NAME = "PhaseIndeterminate";
    private final String MENDELIAN_VIOLATION_FILTER_NAME = "MendelianViolation";
    private final String TRANSMISSION_PROBABILITY_TAG_NAME = "TP";
    private final String SOURCE_NAME = "PhaseByTransmission";

    /**
     * Parse the familial relationship specification, and initialize VCF writer
     */
    public void initialize() {
        String[] pieces = familyStr.split("[\\+\\=]");

        SAMPLE_NAME_MOM = pieces[0];
        SAMPLE_NAME_DAD = pieces[1];
        SAMPLE_NAME_CHILD = pieces[2];

        ArrayList<String> rodNames = new ArrayList<String>();
        rodNames.add(ROD_NAME);

        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        Set<String> vcfSamples = SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        if (vcfSamples.size() != 3) {
            throw new UserException("File to phase by transmission contains more than three samples.  This walker only" +
                                    "accepts VCFs with three samples, so that the meaning of the applied filters is" +
                                    "unambiguous.");
        }

        if (!vcfSamples.contains(SAMPLE_NAME_MOM) || !vcfSamples.contains(SAMPLE_NAME_DAD) || !vcfSamples.contains(SAMPLE_NAME_CHILD)) {
            throw new UserException("One or more of the samples specified in the familyPattern argument is not present" +
                                    "in this file.  Please supply a VCF file that contains only three samples: the" +
                                    "mother, the father, and the child");
        }

        Set<String> samples = new HashSet<String>();
        samples.add(SAMPLE_NAME_MOM);
        samples.add(SAMPLE_NAME_DAD);
        samples.add(SAMPLE_NAME_CHILD);

        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        headerLines.addAll(VCFUtils.getHeaderFields(this.getToolkit()));
        headerLines.add(new VCFFilterHeaderLine(AMBIGUOUS_ALLELE_ORIGIN_FILTER_NAME, "The parental origin of each of the child's allele cannot be determined (ie everyone is heterozygous)"));
        headerLines.add(new VCFFilterHeaderLine(PHASE_INDETERMINATE_FILTER_NAME, "The phase of the child's genotype cannot be determined (ie someone is a no-call)"));
        headerLines.add(new VCFFilterHeaderLine(MENDELIAN_VIOLATION_FILTER_NAME, "No combination of the parents' alleles can yield the child's genotype (ie a possible Mendelian violation)"));
        headerLines.add(new VCFFormatHeaderLine(TRANSMISSION_PROBABILITY_TAG_NAME, 1, VCFHeaderLineType.Float, "Probability that the phase is correct given that the genotypes are correct"));
        vcfWriter.writeHeader(new VCFHeader(headerLines, samples));
    }

    /**
     * Given genotypes for mom, dad, and child, determine the phase for the child by examining the possible transmitted alleles.
     *
     * @param mom    mom's genotype
     * @param dad    dad's genotype
     * @param child  child's genotype
     * @return  phased version of child's genotype
     * @throws  AmbiguousAlleleOriginException if the parentage of the alleles can't be determined (i.e. a triple-het situation),
     * @throws  InsufficientInfoToPhaseGenotypeException if there is insufficient information to determine the phase
     * @throws  MendelianViolationException if no combination of the parental alleles can yield the child's genotype
     */
    private Genotype getPhasedChildGenotype(Genotype mom, Genotype dad, Genotype child) throws AmbiguousAlleleOriginException, InsufficientInfoToPhaseGenotypeException, MendelianViolationException {
        if (mom.isNoCall() || dad.isNoCall() || child.isNoCall()) {
            throw new InsufficientInfoToPhaseGenotypeException();
        }

        Set<Allele> momAlleles = new HashSet<Allele>();
        momAlleles.addAll(mom.getAlleles());

        Set<Allele> dadAlleles = new HashSet<Allele>();
        dadAlleles.addAll(dad.getAlleles());

        double transmissionProb = QualityUtils.qualToProb(mom.getNegLog10PError())*QualityUtils.qualToProb(dad.getNegLog10PError());

        Map<String, Object> attributes = new HashMap<String, Object>();
        attributes.putAll(child.getAttributes());
        attributes.put(TRANSMISSION_PROBABILITY_TAG_NAME, transmissionProb);

        Set<Genotype> possiblePhasedChildGenotypes = new HashSet<Genotype>();

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

        Set<Genotype> ambiguousAlleleOriginGenotypes = new HashSet<Genotype>();

        for (Genotype g1 : possiblePhasedChildGenotypes) {
            for (Genotype g2 : possiblePhasedChildGenotypes) {
                if (!g1.equals(g2)) {
                    if (g1.sameGenotype(g2, true)) {
                        ambiguousAlleleOriginGenotypes.add(g1);
                        ambiguousAlleleOriginGenotypes.add(g2);
                    }
                }
            }
        }

        for (Genotype possiblePhasedChildGenotype : possiblePhasedChildGenotypes) {
            if (child.sameGenotype(possiblePhasedChildGenotype, true)) {
                if (ambiguousAlleleOriginGenotypes.contains(possiblePhasedChildGenotype)) {
                    throw new AmbiguousAlleleOriginException();
                }

                return possiblePhasedChildGenotype;
            }
        }

        throw new MendelianViolationException();
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
            Collection<VariantContext> vcs = tracker.getVariantContexts(ref, ROD_NAME, null, context.getLocation(), true, true);

            for (VariantContext vc : vcs) {
                Genotype mom = vc.getGenotype(SAMPLE_NAME_MOM);
                Genotype dad = vc.getGenotype(SAMPLE_NAME_DAD);
                Genotype child = vc.getGenotype(SAMPLE_NAME_CHILD);

                Set<String> filters = new HashSet<String>();
                filters.addAll(vc.getFilters());

                try {
                    child = getPhasedChildGenotype(mom, dad, child);
                } catch (AmbiguousAlleleOriginException e) {
                    if (!noFilters) {
                        filters.add(AMBIGUOUS_ALLELE_ORIGIN_FILTER_NAME);
                    }
                } catch (InsufficientInfoToPhaseGenotypeException e) {
                    if (!noFilters) {
                        filters.add(PHASE_INDETERMINATE_FILTER_NAME);
                    }
                } catch (MendelianViolationException e) {
                    if (!noFilters) {
                        filters.add(MENDELIAN_VIOLATION_FILTER_NAME);
                    }
                }

                Map<String, Object> attributes = new HashMap<String, Object>();
                attributes.putAll(vc.getAttributes());

                Collection<Genotype> phasedGenotypes = new ArrayList<Genotype>();
                Collection<Allele> alleles = vc.getAlleles();

                phasedGenotypes.add(mom);
                phasedGenotypes.add(dad);
                phasedGenotypes.add(child);

                VariantContext newvc = new VariantContext(SOURCE_NAME, ref.getLocus().getContig(), vc.getStart(), vc.getStart(), alleles, phasedGenotypes, vc.getNegLog10PError(), filters, attributes);

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
