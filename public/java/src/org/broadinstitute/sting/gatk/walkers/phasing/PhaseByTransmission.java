package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.io.File;
import java.io.FileNotFoundException;
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

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Argument(shortName="f", fullName="familySpec", required=true, doc="Patterns for the family structure (usage: mom+dad=child).  Specify several trios by supplying this argument many times and/or a file containing many patterns.")
    public ArrayList<String> familySpecs = null;

    @Output
    protected VCFWriter vcfWriter = null;

    private final String TRANSMISSION_PROBABILITY_TAG_NAME = "TP";
    private final String SOURCE_NAME = "PhaseByTransmission";

    private final Double MENDELIAN_VIOLATION_PRIOR = 1e-8;

    private class Trio {
        private String mother;
        private String father;
        private String child;

        public Trio(String mother, String father, String child) {
            this.mother = mother;
            this.father = father;
            this.child = child;
        }

        public Trio(String familySpec) {
            String[] pieces = familySpec.split("[\\+\\=]");

            this.mother = pieces[0];
            this.father = pieces[1];
            this.child = pieces[2];
        }

        public String getMother() { return mother; }
        public String getFather() { return father; }
        public String getChild() { return child; }
    }

    private ArrayList<Trio> trios = new ArrayList<Trio>();

    public ArrayList<Trio> getFamilySpecsFromCommandLineInput(ArrayList<String> familySpecs) {
        if (familySpecs != null) {
            // Let's first go through the list and see if we were given any files.  We'll add every entry in the file to our
            // spec list set, and treat the entries as if they had been specified on the command line.
            ArrayList<Trio> specs = new ArrayList<Trio>();
            for (String familySpec : familySpecs) {
                File specFile = new File(familySpec);

                try {
                    XReadLines reader = new XReadLines(specFile);

                    List<String> lines = reader.readLines();
                    for (String line : lines) {
                        specs.add(new Trio(line));
                    }
                } catch (FileNotFoundException e) {
                    specs.add(new Trio(familySpec)); // not a file, so must be a family spec
                }
            }

            return specs;
        }

        return new ArrayList<Trio>();
    }

    /**
     * Parse the familial relationship specification, and initialize VCF writer
     */
    public void initialize() {
        trios = getFamilySpecsFromCommandLineInput(familySpecs);

        ArrayList<String> rodNames = new ArrayList<String>();
        rodNames.add(variantCollection.variants.getName());

        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        Set<String> vcfSamples = SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        headerLines.addAll(VCFUtils.getHeaderFields(this.getToolkit()));
        headerLines.add(new VCFFormatHeaderLine(TRANSMISSION_PROBABILITY_TAG_NAME, 1, VCFHeaderLineType.Float, "Probability that the phase is correct given that the genotypes are correct"));
        headerLines.add(new VCFHeaderLine("source", SOURCE_NAME));
        vcfWriter.writeHeader(new VCFHeader(headerLines, vcfSamples));
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

    private ArrayList<Genotype> phaseTrioGenotypes(Allele ref, Allele alt, Genotype mother, Genotype father, Genotype child) {
        ArrayList<Genotype> finalGenotypes = new ArrayList<Genotype>();
        finalGenotypes.add(mother);
        finalGenotypes.add(father);
        finalGenotypes.add(child);

        if (mother.isCalled() && father.isCalled() && child.isCalled()) {
            ArrayList<Genotype> possibleMotherGenotypes = createAllThreeGenotypes(ref, alt, mother);
            ArrayList<Genotype> possibleFatherGenotypes = createAllThreeGenotypes(ref, alt, father);
            ArrayList<Genotype> possibleChildGenotypes = createAllThreeGenotypes(ref, alt, child);

            double bestConfigurationLikelihood = 0.0;
            double bestPrior = 0.0;
            Genotype bestMotherGenotype = mother;
            Genotype bestFatherGenotype = father;
            Genotype bestChildGenotype = child;

            double norm = 0.0;

            for (Genotype motherGenotype : possibleMotherGenotypes) {
                for (Genotype fatherGenotype : possibleFatherGenotypes) {
                    for (Genotype childGenotype : possibleChildGenotypes) {
                        double prior = isMendelianViolation(ref, alt, motherGenotype, fatherGenotype, childGenotype) ? MENDELIAN_VIOLATION_PRIOR : 1.0 - 12*MENDELIAN_VIOLATION_PRIOR;
                        double configurationLikelihood = computeTransmissionLikelihoodOfGenotypeConfiguration(motherGenotype, fatherGenotype, childGenotype);
                        norm += prior*configurationLikelihood;

                        if (prior*configurationLikelihood > bestPrior*bestConfigurationLikelihood) {
                            bestConfigurationLikelihood = configurationLikelihood;
                            bestPrior = prior;
                            bestMotherGenotype = motherGenotype;
                            bestFatherGenotype = fatherGenotype;
                            bestChildGenotype = childGenotype;
                        }
                    }
                }
            }

            if (!(bestMotherGenotype.isHet() && bestFatherGenotype.isHet() && bestChildGenotype.isHet())) {
                Map<String, Object> attributes = new HashMap<String, Object>();
                attributes.putAll(bestChildGenotype.getAttributes());
                attributes.put(TRANSMISSION_PROBABILITY_TAG_NAME, bestPrior*bestConfigurationLikelihood / norm);
                bestChildGenotype = Genotype.modifyAttributes(bestChildGenotype, attributes);

                finalGenotypes = getPhasedGenotypes(bestMotherGenotype, bestFatherGenotype, bestChildGenotype);
            }
        }

        return finalGenotypes;
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
            VariantContext vc = tracker.getFirstValue(variantCollection.variants, context.getLocation());

            Map<String, Genotype> genotypeMap = vc.getGenotypes();

            for (Trio trio : trios) {
                Genotype mother = vc.getGenotype(trio.getMother());
                Genotype father = vc.getGenotype(trio.getFather());
                Genotype child = vc.getGenotype(trio.getChild());

                ArrayList<Genotype> trioGenotypes = phaseTrioGenotypes(vc.getReference(), vc.getAltAlleleWithHighestAlleleCount(), mother, father, child);

                Genotype phasedMother = trioGenotypes.get(0);
                Genotype phasedFather = trioGenotypes.get(1);
                Genotype phasedChild = trioGenotypes.get(2);

                genotypeMap.put(phasedMother.getSampleName(), phasedMother);
                genotypeMap.put(phasedFather.getSampleName(), phasedFather);
                genotypeMap.put(phasedChild.getSampleName(), phasedChild);
            }

            VariantContext newvc = VariantContext.modifyGenotypes(vc, genotypeMap);

            vcfWriter.add(newvc);
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
