package org.broadinstitute.sting.playground.gatk.walkers.replication_validation;

import com.google.common.collect.ArrayListMultimap;
import com.sun.xml.internal.ws.client.BindingProviderProperties;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.MathUtils;

import java.io.PrintStream;
import org.apache.commons.math.util.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.*;
import java.util.regex.Pattern;

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
     * @param data - array of bases from the reference sample
     * @param refBase - base from the reference sequence for this site
     * @param phredScaledPrior - phred scaled expected site quality (prior)
     * @return an array of log10 probabilities of site qualities ranging from Q1-Q40.
     */
    private double[] buildErrorModel (byte[] data, byte refBase, byte phredScaledPrior) {
        double [] model = new double[maxQualityScore+1];
        int coverage = data.length;
        int mismatches = getNumberOfMismatches(data, refBase);
        int matches = coverage - mismatches;

        for (byte q=0; q<=maxQualityScore; q++) {
            double probMismatch = MathUtils.phredScaleToProbability(q);
            model[q] = MathUtils.phredScaleToLog10Probability(phredScaledPrior) +
                       org.apache.commons.math.util.MathUtils.binomialCoefficientLog(coverage, mismatches) +
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
    private int getNumberOfMismatches(byte[] data, byte refBase) {
        int mismatches = 0;
        for (byte seqBase : data) {
            if (seqBase != refBase)
                mismatches++;
        }
        return mismatches;
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
        VariantContext referenceContext = tracker.getVariantContext(ref, REFERENCE_ROD_NAME, context.getLocation());
        byte trueReferenceBase;

        // Site is not a variant, take from the reference
        if (referenceContext == null) {
            trueReferenceBase = ref.getBase();
        }

        // Site has a VCF entry -- is variant
        else {
            // Site is filtered, don't mess with it if option is set
            if (referenceContext.isFiltered() && EXCLUDE_FILTERED_REFERENCE_SITES) {
                return 0;
            }

            // Site is HomVar -- easy, take the only alt allele.
            else if (referenceContext.getGenotype(referenceSampleName).isHomVar()) {
                trueReferenceBase = referenceContext.getAlternateAllele(0).getBases()[0];
            }

            // TODO: treat HET cases specially?
            else {
                return 0;
            }
        }



        ReadBackedPileup contextPileup = context.getBasePileup();
        Set<String> lanesInLocus = getLaneIDs(contextPileup.getReadGroups());
        for (String laneID : lanesInLocus) {
            // make a pileup for this lane
            ReadBackedPileup lanePileup = contextPileup.getPileupForLane(laneID);

            // make a pileup without the reference sample
            HashSet<String> samples = new HashSet<String>(lanePileup.getSampleNames());
            samples.remove(referenceSampleName);

            ReadBackedPileup poolPileup = lanePileup.getPileupForSampleNames(samples);

            //System.out.println(ref.getLocus().getContig() + ":" + ref.getLocus().getStart());

            // create a reference pileup to build error model
            ReadBackedPileup referenceSamplePileup = lanePileup.getPileupForSampleName(referenceSampleName);

            // We can only analyze loci where we have reference sample coverage. And it only makes sense to build
            // an error model if the site has reads in the pileup (other than the reference sample)
            if (referenceSamplePileup != null && !referenceSamplePileup.isEmpty() && !poolPileup.isEmpty()) {

                // Build error model
                double [] errorModel = buildErrorModel(referenceSamplePileup.getBases(), trueReferenceBase, defaultPrior);

                // Debug error model
                System.out.println(laneID + ":");
                for (double v : errorModel)
                    System.out.print(v + ", ");
                System.out.println();

                // todo: call each pool for this site
                // todo: merge pools
                // todo: decide whether or not it's a variant
            }
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

