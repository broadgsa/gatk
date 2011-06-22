package org.broadinstitute.sting.playground.gatk.walkers.replication_validation;

import com.google.common.collect.ArrayListMultimap;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.MathUtils;

import java.io.PrintStream;
import org.apache.commons.math.util.*;
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
    @Argument(shortName="refSample", fullName="reference_sample_name", doc="reference sample name", required=true)
    String referenceSampleName = "NA12878";

    @Argument(shortName="nsamples", fullName="number_of_samples", doc="Number of samples in the dataset (not counting the reference sample).", required=true)
    int nSamples = -1;

    @Argument(shortName="nchr", fullName="number_of_chromosomes", doc="Number of chromosomes per sample (in case you're not dealing with diploids). Default: 2.", required=false)
    int nChromosomes = 2;

    @Argument(shortName="maxac", fullName="max_allele_count", doc="Max number of alleles expected in a site. Default: 2 * number of samples. Smaller numbers process faster.", required=false)
    int overrideMaxAlleleCount = -1;

    @Argument(shortName="maxqs", fullName="max_quality_score", doc="Max quality score to consider. Default: Q40. Smaller numbers process faster.", required=false)
    int maxQualityScore= 40;

    @Argument(shortName="prior", fullName="site_quality_prior", doc="Phred-Scaled prior quality of the site. Default: Q20.", required=false)
    byte defaultPrior= 20;

    @Output(doc="Write output to this file instead of STDOUT")
    PrintStream out;

    int maxAlleleCount;
    List<LaneData> lanes;

    private class LaneData {
        String id;
        List <String> sampleNames;
        String referenceSample;

        public String getId() {
            return id;
        }

        public void setId(String id) {
            this.id = id;
        }

        public List<String> getSampleNames() {
            return sampleNames;
        }

        public void setSampleNames(List<String> sampleNames) {
            this.sampleNames = sampleNames;
        }

        public String getReferenceSample() {
            return referenceSample;
        }

        public void setReferenceSample(String referenceSample) {
            this.referenceSample = referenceSample;
        }
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

    // buildErrorModel helper functions
    private int getNumberOfMismatches(byte[] data, byte refBase) {
        int mismatches = 0;
        for (byte seqBase : data) {
            if (seqBase != refBase)
                mismatches++;
        }
        return mismatches;
    }

//    private String getLaneId (SAMReadGroupRecord rg) {
//        String [] idParts = rg.getPlatformUnit().split("\\.");
//        return idParts[0] + idParts[1];  // 0 is run_id, 1 is lane_number.
//    }
//
//    private List<LaneData> createLanesList () {
//        List<SAMReadGroupRecord> rgs = getToolkit().getSAMFileHeader().getReadGroups();
//        for(SAMReadGroupRecord rg : rgs) {
//            String laneId = getLaneId(rg);
//        }
//
//        return new LinkedList<LaneData>();
//    }
//
//    private int getNumberOfSamples() {
//        int nSamples = 0;
//        List<Set<String>> readers = getToolkit().getSamplesByReaders();
//        for(Set<String> reader : readers) {
//            nSamples += reader.size();
//            if (reader.contains(referenceSampleName))
//                nSamples--;
//        }
//        return nSamples;
//    }


    public void initialize() {
        if (overrideMaxAlleleCount > 0)
            maxAlleleCount = overrideMaxAlleleCount;
        else
            maxAlleleCount = nSamples*nChromosomes;
    }



    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
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
                double [] errorModel = buildErrorModel(poolPileup.getBases(), referenceSamplePileup.getBases()[0], defaultPrior);

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

    /**
     * GATK Engine creates readgroups of the form XXX.Y.Z
     * XXX.Y is the unique lane identifier.
     *     Z is the id of the sample to make the read group id unique
     * This function returns a map that groups together all samples that were run in the same lane:
     *     XXX.Y -> {Z1, Z2, Z3, ... Zn}
     * @param readGroups A collection of read group strings (obtained from the alignment context pileup)
     * @return A map of lanes and samples, as described above
     */
//    private HashMap<String, ArrayList<String>> decodeReadIDs(Collection<String> readGroups) {
//        HashMap<String, ArrayList<String>> result = new HashMap<String, ArrayList<String>>();
//        for (String rgid : readGroups) {
//            String [] parsedId = rgid.split("\\.");
//            String lane = parsedId[0] + "." + parsedId[1];
//            String id = parsedId[2];
//            ArrayList<String> idList = (result.get(lane) == null) ? (new ArrayList<String>()) : result.get(lane);
//            idList.add(id);
//            result.put(lane, idList);
//        }
//        return result;
//    }

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


    public Long reduceInit() {
        return 0l;
    }

    public Long reduce(Integer value, Long sum) {
        return value + sum;
    }

    /**
     * Reduces two subtrees together.  In this case, the implementation of the tree reduce
     * is exactly the same as the implementation of the single reduce.
     */
    public Long treeReduce(Long lhs, Long rhs) {
        return lhs + rhs;
    }

    public void onTraversalDone( Long c ) {
        out.println(c);
    }
}

