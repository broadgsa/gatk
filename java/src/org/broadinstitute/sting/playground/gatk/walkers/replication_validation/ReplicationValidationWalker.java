package org.broadinstitute.sting.playground.gatk.walkers.replication_validation;

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.MathUtils;

import java.io.PrintStream;

import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.*;

/**
 * Implementation of the replication and validation framework with reference based error model
 * for pooled sequencing.
 *
 * The input should be a BAM file with pooled sequencing data where each pool is represented by
 * samples with the same barcode.
 *
 * A reference sample name must be provided and it must be barcoded uniquely.
 */
public class ReplicationValidationWalker extends LocusWalker<Integer, Long> implements TreeReducible<Long> {


    @Argument(shortName="refsample", fullName="reference_sample_name", doc="Reference sample name.", required=true)
    String referenceSampleName;

    @Argument(shortName="nsamples", fullName="number_of_samples", doc="Number of samples in the dataset (not counting the reference sample).", required=true)
    int nSamples = -1;

    @Argument(shortName="nchr", fullName="number_of_chromosomes", doc="Number of chromosomes per sample (in case you're not dealing with diploids). Default: 2.", required=false)
    int nChromosomes = 2;

    @Argument(shortName="maxac", fullName="max_allele_count", doc="Max number of alleles expected in a site. Smaller numbers process faster. Default: 2 * number of samples. ", required=false)
    int overrideMaxAlleleCount = -1;

    @Argument(shortName="maxqs", fullName="max_quality_score", doc="Max quality score to consider. Smaller numbers process faster. Default: Q40.", required=false)
    int maxQualityScore= 40;

    @Argument(shortName="prior", fullName="site_quality_prior", doc="Phred-Scaled prior quality of the site. Default: Q20.", required=false)
    byte defaultPrior= 20;

    @Argument(shortName="ef", fullName="exclude_filtered_reference_sites", doc="Don't include in the analysis sites where the reference sample VCF is filtered. Default: false.", required=false)
    boolean EXCLUDE_FILTERED_REFERENCE_SITES = false;


    @Output(doc="Write output to this file instead of STDOUT")
    PrintStream out;

    int maxAlleleCount;
    boolean USE_TRUTH_ROD;

    final String REFERENCE_ROD_NAME = "reference";
    final String TRUTH_ROD_NAME = "truth";


    /**
     * GATK Engine creates readgroups of the form XXX.Y.Z
     * XXX.Y is the unique lane identifier.
     *     Z is the id of the sample to make the read group id unique
     * This function returns a list of unique lane identifiers.
     * @param readGroups readGroups A collection of read group strings (obtained from the alignment context pileup)
     * @return a collection of lane ids.
     */
    private Set<String> getLaneIDs(Collection<String> readGroups) {
        HashSet<String> result = new HashSet<String>();
        for (String rgid : readGroups) {
            String [] parsedId = rgid.split("\\.");
            result.add(parsedId[0] + "." + parsedId[1]);
        }
        return result;
    }

    /**
     * Calculates he probability of the data (reference sample reads) given the phred scaled site quality score.
     * @param referenceSamplePileup reference sample pileup
     * @param refBases base from the reference sequence for this site
     * @param phredScaledPrior phred scaled expected site quality (prior)
     * @return an array of log10 probabilities of site qualities ranging from Q1-Q40.
     */
    private double[] buildErrorModel (ReadBackedPileup referenceSamplePileup, Collection<Byte> refBases, byte phredScaledPrior) {
        double [] model = new double[maxQualityScore+1];
        byte [] data = referenceSamplePileup.getBases();

        int coverage = data.length;
        int mismatches = getNumberOfMismatches(data, refBases);
        int matches = coverage - mismatches;

        for (byte q=0; q<=maxQualityScore; q++) {
            double probMismatch = MathUtils.phredScaleToProbability(q);
            model[q] = MathUtils.phredScaleToLog10Probability(phredScaledPrior) +
                       MathUtils.log10BinomialCoefficient(coverage, mismatches) +
                       mismatches * Math.log10(probMismatch) +
                       matches * Math.log10(1-probMismatch);
        }
        return model;
    }

    /**
     * Returns the number of mismatching bases in a pileup
     * @param data the bases of a pileup
     * @param refBase the reference sample base to compare to
     * @return number of bases in data that are different from refBase
     */
    private int getNumberOfMismatches (byte[] data, Collection<Byte> refBase) {
        int mismatches = 0;
        for (byte seqBase : data) {
            if (!refBase.contains(seqBase))
                mismatches++;
        }
        return mismatches;
    }

    /**
     * Returns the true bases for the reference sample in this locus. Homozygous loci will return one base
     * but heterozygous will return two bases (hence why it returns a collection).
     *
     * @param referenceSampleContext the variant context from the reference sample ROD track
     * @param ref the reference sequence context
     * @return the true bases for the reference sample.
     */
    private Collection<Byte> getTrueBases(VariantContext referenceSampleContext, ReferenceContext ref) {

        ArrayList<Byte> trueReferenceBase = new ArrayList<Byte>();

        // Site is not a variant, take from the reference
        if (referenceSampleContext == null) {
            trueReferenceBase.add(ref.getBase());
        }

        else if (referenceSampleContext.isIndel()) {
            return null; // TODO: add special treatment for extended events. For Now just skip these altogether.
        }

        // Site has a VCF entry -- is variant
        else {
            // Site is filtered, don't mess with it if option is set
            if (referenceSampleContext.isFiltered() && EXCLUDE_FILTERED_REFERENCE_SITES) {
                return null;
            }

            Genotype referenceGenotype = referenceSampleContext.getGenotype(referenceSampleName);
            List<Allele> referenceAlleles = referenceGenotype.getAlleles();
            for (Allele allele : referenceAlleles) {
                byte [] bases = allele.getBases();
                for (byte b : bases) {
                    trueReferenceBase.add(b);
                }
            }
        }
        return trueReferenceBase;
    }

    public void initialize() {

        // Set the max allele count (defines the size of the error model array)
        maxAlleleCount = (overrideMaxAlleleCount > 0) ? overrideMaxAlleleCount : nSamples*nChromosomes;

        // Look for the reference ROD and the optional truth ROD. If truth is provided, set the truth "test" mode ON.
        List<ReferenceOrderedDataSource> rods = getToolkit().getRodDataSources();
        if (rods.size() < 1) {
            throw new IllegalArgumentException("You must provide a reference ROD.");
        }
        boolean foundReferenceROD = false;
        boolean foundTruthROD = false;
        for (ReferenceOrderedDataSource rod : rods) {
            if (rod.getName().equals(REFERENCE_ROD_NAME)) {
                foundReferenceROD = true;
            }
            if (rod.getName().equals(TRUTH_ROD_NAME)) {
                foundReferenceROD = true;
            }
        }
        if (!foundReferenceROD) {
            throw new IllegalArgumentException("You haven't provided a reference ROD. Note that the reference ROD must be labeled " + REFERENCE_ROD_NAME + ".");
        }
        if (rods.size() > 1 && !foundTruthROD) {
            throw new IllegalArgumentException("You haven't provided a truth ROD. Note that the reference ROD must be labeled " + TRUTH_ROD_NAME + ".");
        }
        USE_TRUTH_ROD = foundTruthROD;
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        // Get reference base from VCF or Reference
        VariantContext referenceSampleContext = tracker.getVariantContext(ref, REFERENCE_ROD_NAME, context.getLocation());
        Collection<Byte> trueReferenceBases = getTrueBases(referenceSampleContext, ref);

        // If there is no true reference base in this locus, skip it.
        if (trueReferenceBases == null)
            return 0;

        ReadBackedPileup contextPileup = context.getBasePileup();
        Set<String> lanesInLocus = getLaneIDs(contextPileup.getReadGroups());
        for (String laneID : lanesInLocus) {
            // make a pileup for this lane
            ReadBackedPileup lanePileup = contextPileup.getPileupForLane(laneID);
            Collection<String> samplesInLane = lanePileup.getSampleNames();

            // we can only analyze loci that have reads for the reference sample
            if (samplesInLane.contains(referenceSampleName)) {

                // build reference sample pileup
                ReadBackedPileup referenceSamplePileup = lanePileup.getPileupForSampleName(referenceSampleName);

                // Build error model
                double [] errorModel = buildErrorModel(referenceSamplePileup, trueReferenceBases, defaultPrior);

                // iterate over all samples (pools) in this lane except the reference
                samplesInLane.remove(referenceSampleName);
                for (String pool : samplesInLane) {
                    ReadBackedPileup poolPileup = lanePileup.getPileupForSampleName(pool);

                    // Debug error model
                    if (referenceSamplePileup.getBases().length > 50) {
                        System.out.println("\n" + laneID + " - " + pool + ": " + referenceSamplePileup.getBases().length);
                        for (double v : errorModel)
                            System.out.print(v + ", ");
                        System.out.println();
                    }
                }
            }
            // todo: call each pool for this site
            // todo: merge pools
            // todo: decide whether or not it's a variant
        }
        return 1;
    }

    public Long reduceInit() {
        return 0l;
    }

    public Long reduce(Integer value, Long sum) {
        return value + sum;
    }

    public Long treeReduce(Long lhs, Long rhs) {
        return lhs + rhs;
    }

    public void onTraversalDone( Long c ) {
        out.println(c);
    }
}

