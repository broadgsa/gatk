package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.util.*;

/**
 * Phases a trio VCF (child phased by transmission, implied phase carried over to parents).  Given genotypes for a trio,
 * this walker modifies the genotypes (if necessary) to reflect the most likely configuration given the genotype
 * likelihoods and inheritance constraints, phases child by transmission and carries over implied phase to the parents
 * (their alleles in their genotypes are ordered as transmitted|untransmitted).  Computes probability that the
 * determined phase is correct given that the genotype configuration is correct (useful if you want to use this to
 * compare phasing accuracy, but want to break that comparison down by phasing confidence in the truth set).  Optionally
 * filters out sites where the phasing is indeterminate (site has no-calls), ambiguous (everyone is heterozygous), or
 * the genotypes exhibit a Mendelian violation.  This walker assumes there are only three samples in the VCF file to
 * begin.
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

    private final String ROD_NAME = "variant";
    private final String AMBIGUOUS_ALLELE_ORIGIN_FILTER_NAME = "AmbiguousAlleleOrigin";
    private final String INSUFFICIENT_DATA_FILTER_NAME = "InsufficientInformation";
    private final String MENDELIAN_VIOLATION_FILTER_NAME = "MendelianViolation";
    private final String TRANSMISSION_PROBABILITY_TAG_NAME = "TP";
    private final String SOURCE_NAME = "PhaseByTransmission";

    private final Double MENDELIAN_VIOLATION_PRIOR = 1e-8;

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

        Set<String> samples = new TreeSet<String>();
        samples.add(SAMPLE_NAME_MOM);
        samples.add(SAMPLE_NAME_DAD);
        samples.add(SAMPLE_NAME_CHILD);

        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        headerLines.addAll(VCFUtils.getHeaderFields(this.getToolkit()));

        if (!noFilters) {
            headerLines.add(new VCFFilterHeaderLine(AMBIGUOUS_ALLELE_ORIGIN_FILTER_NAME, "The parental origin of each of the child's allele cannot be determined (ie everyone is heterozygous)"));
            headerLines.add(new VCFFilterHeaderLine(INSUFFICIENT_DATA_FILTER_NAME, "The phase of the child's genotype cannot be determined (ie someone is a no-call)"));
            headerLines.add(new VCFFilterHeaderLine(MENDELIAN_VIOLATION_FILTER_NAME, "No combination of the parents' alleles can yield the child's genotype (ie a possible Mendelian violation)"));
        }

        headerLines.add(new VCFInfoHeaderLine(TRANSMISSION_PROBABILITY_TAG_NAME, 1, VCFHeaderLineType.Float, "Probability that the phase is correct given that the genotypes are correct"));
        vcfWriter.writeHeader(new VCFHeader(headerLines, samples));
    }

    private double computeTransmissionLikelihoodOfGenotypeConfiguration(Genotype mom, Genotype dad, Genotype child) {
        double[] momLikelihoods = MathUtils.normalizeFromLog10(mom.getLikelihoods().getAsVector());
        double[] dadLikelihoods = MathUtils.normalizeFromLog10(dad.getLikelihoods().getAsVector());
        double[] childLikelihoods = MathUtils.normalizeFromLog10(child.getLikelihoods().getAsVector());

        int momIndex = mom.getType().ordinal() - 1;
        int dadIndex = dad.getType().ordinal() - 1;
        int childIndex = child.getType().ordinal() - 1;

        return momLikelihoods[momIndex]*dadLikelihoods[dadIndex]*childLikelihoods[childIndex];
    }

    private ArrayList<Genotype> createAllThreeGenotypes(Allele refAllele, Allele altAllele, Genotype g) {
        List<Allele> homRefAlleles = new ArrayList<Allele>();
        homRefAlleles.add(refAllele);
        homRefAlleles.add(refAllele);
        Genotype homRef = new Genotype(g.getSampleName(), homRefAlleles, g.getNegLog10PError(), null, g.getAttributes(), false);

        List<Allele> hetAlleles = new ArrayList<Allele>();
        hetAlleles.add(refAllele);
        hetAlleles.add(altAllele);
        Genotype het = new Genotype(g.getSampleName(), hetAlleles, g.getNegLog10PError(), null, g.getAttributes(), false);

        List<Allele> homVarAlleles = new ArrayList<Allele>();
        homVarAlleles.add(altAllele);
        homVarAlleles.add(altAllele);
        Genotype homVar = new Genotype(g.getSampleName(), homVarAlleles, g.getNegLog10PError(), null, g.getAttributes(), false);

        ArrayList<Genotype> genotypes = new ArrayList<Genotype>();
        genotypes.add(homRef);
        genotypes.add(het);
        genotypes.add(homVar);

        return genotypes;
    }

    private int getNumberOfMatchingAlleles(Allele alleleToMatch, Genotype g) {
        List<Allele> alleles = g.getAlleles();
        int matchingAlleles = 0;

        for (Allele a : alleles) {
            if (!alleleToMatch.equals(a)) {
                matchingAlleles++;
            }
        }

        return matchingAlleles;
    }

    private boolean isMendelianViolation(Allele refAllele, Allele altAllele, Genotype mom, Genotype dad, Genotype child) {
        int numMomRefAlleles = getNumberOfMatchingAlleles(refAllele, mom) > 0 ? 1 : 0;
        int numMomAltAlleles = getNumberOfMatchingAlleles(altAllele, mom) > 0 ? 1 : 0;

        int numDadRefAlleles = getNumberOfMatchingAlleles(refAllele, dad) > 0 ? 1 : 0;
        int numDadAltAlleles = getNumberOfMatchingAlleles(altAllele, dad) > 0 ? 1 : 0;

        int numChildRefAlleles = getNumberOfMatchingAlleles(refAllele, child);
        int numChildAltAlleles = getNumberOfMatchingAlleles(altAllele, child);

        return (numMomRefAlleles + numDadRefAlleles < numChildRefAlleles || numMomAltAlleles + numDadAltAlleles < numChildAltAlleles);
    }

    private ArrayList<Genotype> getPhasedGenotypes(Genotype mom, Genotype dad, Genotype child) {
        Set<Genotype> possiblePhasedChildGenotypes = new HashSet<Genotype>();

        for (Allele momAllele : mom.getAlleles()) {
            for (Allele dadAllele : dad.getAlleles()) {
                ArrayList<Allele> possiblePhasedChildAlleles = new ArrayList<Allele>();
                possiblePhasedChildAlleles.add(momAllele);
                possiblePhasedChildAlleles.add(dadAllele);

                Genotype possiblePhasedChildGenotype = new Genotype(child.getSampleName(), possiblePhasedChildAlleles, child.getNegLog10PError(), child.getFilters(), child.getAttributes(), true);

                possiblePhasedChildGenotypes.add(possiblePhasedChildGenotype);
            }
        }

        ArrayList<Genotype> finalGenotypes = new ArrayList<Genotype>();

        for (Genotype phasedChildGenotype : possiblePhasedChildGenotypes) {
            if (child.sameGenotype(phasedChildGenotype, true)) {
                Allele momTransmittedAllele = phasedChildGenotype.getAllele(0);
                Allele momUntransmittedAllele = mom.getAllele(0) != momTransmittedAllele ? mom.getAllele(0) : mom.getAllele(1);

                ArrayList<Allele> phasedMomAlleles = new ArrayList<Allele>();
                phasedMomAlleles.add(momTransmittedAllele);
                phasedMomAlleles.add(momUntransmittedAllele);

                Genotype phasedMomGenotype = new Genotype(mom.getSampleName(), phasedMomAlleles, mom.getNegLog10PError(), mom.getFilters(), mom.getAttributes(), true);

                Allele dadTransmittedAllele = phasedChildGenotype.getAllele(1);
                Allele dadUntransmittedAllele = dad.getAllele(0) != dadTransmittedAllele ? dad.getAllele(0) : dad.getAllele(1);

                ArrayList<Allele> phasedDadAlleles = new ArrayList<Allele>();
                phasedDadAlleles.add(dadTransmittedAllele);
                phasedDadAlleles.add(dadUntransmittedAllele);

                Genotype phasedDadGenotype = new Genotype(dad.getSampleName(), phasedDadAlleles, dad.getNegLog10PError(), dad.getFilters(), dad.getAttributes(), true);

                finalGenotypes.add(phasedMomGenotype);
                finalGenotypes.add(phasedDadGenotype);
                finalGenotypes.add(phasedChildGenotype);

                return finalGenotypes;
            }
        }

        finalGenotypes.add(mom);
        finalGenotypes.add(dad);
        finalGenotypes.add(child);

        return finalGenotypes;
    }

    private VariantContext phaseTrioGenotypes(VariantContext vc) {
        Genotype mom = vc.getGenotype(SAMPLE_NAME_MOM);
        Genotype dad = vc.getGenotype(SAMPLE_NAME_DAD);
        Genotype child = vc.getGenotype(SAMPLE_NAME_CHILD);

        Set<String> filters = new HashSet<String>();
        filters.addAll(vc.getFilters());

        Map<String, Object> attributes = new HashMap<String, Object>();
        attributes.putAll(vc.getAttributes());
        attributes.put(TRANSMISSION_PROBABILITY_TAG_NAME, 0.0);

        ArrayList<Genotype> finalGenotypes = new ArrayList<Genotype>();
        finalGenotypes.add(mom);
        finalGenotypes.add(dad);
        finalGenotypes.add(child);

        if (!mom.isCalled() || !dad.isCalled() || !child.isCalled()) {
            filters.add(INSUFFICIENT_DATA_FILTER_NAME);
        } else {
            ArrayList<Genotype> possibleMomGenotypes = createAllThreeGenotypes(vc.getReference(), vc.getAlternateAllele(0), mom);
            ArrayList<Genotype> possibleDadGenotypes = createAllThreeGenotypes(vc.getReference(), vc.getAlternateAllele(0), dad);
            ArrayList<Genotype> possibleChildGenotypes = createAllThreeGenotypes(vc.getReference(), vc.getAlternateAllele(0), child);

            double bestConfigurationLikelihood = 0.0;
            double bestPrior = 0.0;
            Genotype bestMomGenotype = mom;
            Genotype bestDadGenotype = dad;
            Genotype bestChildGenotype = child;

            double norm = 0.0;

            for (Genotype momGenotype : possibleMomGenotypes) {
                for (Genotype dadGenotype : possibleDadGenotypes) {
                    for (Genotype childGenotype : possibleChildGenotypes) {
                        double prior = isMendelianViolation(vc.getReference(), vc.getAlternateAllele(0), momGenotype, dadGenotype, childGenotype) ? MENDELIAN_VIOLATION_PRIOR : 1.0 - 12*MENDELIAN_VIOLATION_PRIOR;
                        double configurationLikelihood = computeTransmissionLikelihoodOfGenotypeConfiguration(momGenotype, dadGenotype, childGenotype);
                        norm += prior*configurationLikelihood;

                        if (prior*configurationLikelihood > bestPrior*bestConfigurationLikelihood) {
                            bestConfigurationLikelihood = configurationLikelihood;
                            bestPrior = prior;
                            bestMomGenotype = momGenotype;
                            bestDadGenotype = dadGenotype;
                            bestChildGenotype = childGenotype;
                        }
                    }
                }
            }

            if (isMendelianViolation(vc.getReference(), vc.getAlternateAllele(0), bestMomGenotype, bestDadGenotype, bestChildGenotype)) {
                filters.add(MENDELIAN_VIOLATION_FILTER_NAME);
            } else if (bestMomGenotype.isHet() && bestDadGenotype.isHet() && bestChildGenotype.isHet()) {
                filters.add(AMBIGUOUS_ALLELE_ORIGIN_FILTER_NAME);
            } else {
                finalGenotypes = getPhasedGenotypes(bestMomGenotype, bestDadGenotype, bestChildGenotype);

                attributes.put(TRANSMISSION_PROBABILITY_TAG_NAME, bestPrior*bestConfigurationLikelihood / norm);
            }
        }

        return new VariantContext(SOURCE_NAME, vc.getChr(), vc.getStart(), vc.getStart(), vc.getAlleles(), finalGenotypes, vc.getNegLog10PError(), noFilters ? vc.getFilters() : filters, attributes);
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
                vcfWriter.add(phaseTrioGenotypes(vc), ref.getBase());
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
