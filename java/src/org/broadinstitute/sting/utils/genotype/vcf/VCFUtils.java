package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.Genotype;

import java.util.*;

/**
 * A set of static utility methods for common operations on VCF files/records.
 */
public class VCFUtils {
    /**
     * Constructor access disallowed...static utility methods only!
     */
    private VCFUtils() { }


    /**
     * Gets the sample names from all VCF rods input by the user and uniquifies them if there is overlap
     * (e.g. sampleX.1, sampleX.2, ...)
     * When finished, samples contains the uniquified sample names and rodNamesToSampleNames contains a mapping
     * from rod/sample pairs to the new uniquified names
     * 
     * @param toolkit    GATK engine
     * @param samples    set to store the sample names
     * @param rodNamesToSampleNames mapping of rod/sample pairs to new uniquified sample names
     */
    public static void getUniquifiedSamplesFromRods(GenomeAnalysisEngine toolkit, Set<String> samples, Map<Pair<String, String>, String> rodNamesToSampleNames) {

        // keep a map of sample name to next available uniquified index
        HashMap<String, Integer> sampleOverlapMap = new HashMap<String, Integer>();

        // iterate to get all of the sample names
        List<ReferenceOrderedDataSource> dataSources = toolkit.getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            ReferenceOrderedData rod = source.getReferenceOrderedData();
            if ( rod.getType().equals(RodVCF.class) ) {
                VCFReader reader = new VCFReader(rod.getFile());
                Set<String> vcfSamples = reader.getHeader().getGenotypeSamples();
                for ( String sample : vcfSamples )
                    addUniqueSample(samples, sampleOverlapMap, rodNamesToSampleNames, sample, rod.getName());
                reader.close();
            }
        }
    }

    private static void addUniqueSample(Set<String> samples, Map<String, Integer> sampleOverlapMap, Map<Pair<String, String>, String> rodNamesToSampleNames, String newSample, String rodName) {

        // if it's already a non-unique sample name, give it a unique suffix and increment the value
        Integer uniqueIndex = sampleOverlapMap.get(newSample);
        if ( uniqueIndex != null ) {
            String uniqueName = newSample + "." + uniqueIndex;
            samples.add(uniqueName);
            rodNamesToSampleNames.put(new Pair<String, String>(rodName, newSample), uniqueName);
            sampleOverlapMap.put(newSample, uniqueIndex + 1);
        }

        // if this is the second occurrence of the sample name, uniquify both of them
        else if ( samples.contains(newSample) ) {
            samples.remove(newSample);
            String uniqueName1 = newSample + "." + 1;
            samples.add(uniqueName1);
            for ( java.util.Map.Entry<Pair<String, String>, String> entry : rodNamesToSampleNames.entrySet() ) {
                if ( entry.getValue().equals(newSample) ) {
                    entry.setValue(uniqueName1);
                    break;
                }
            }

            String uniqueName2 = newSample + "." + 2;
            samples.add(uniqueName2);
            rodNamesToSampleNames.put(new Pair<String, String>(rodName, newSample), uniqueName2);

            sampleOverlapMap.put(newSample, 3);
        }

        // otherwise, just add it to the list of samples
        else {
            samples.add(newSample);
            rodNamesToSampleNames.put(new Pair<String, String>(rodName, newSample), newSample);
        }
    }

    /**
     * Merges various vcf records into a single one using the mapping from rodNamesToSampleNames to get unique sample names
     *
     * @param rods                   the vcf rods
     * @param rodNamesToSampleNames  mapping of rod/sample pairs to new uniquified sample names
     * @return the new merged vcf record
     */
    public static VCFRecord mergeRecords(List<RodVCF> rods, Map<Pair<String, String>, String> rodNamesToSampleNames) {

        VCFParameters params = new VCFParameters();
        params.addFormatItem("GT");

        // keep track of the locus specific data so we can merge them intelligently
        int totalReadDepth = 0;
        double maxConfidence = 0.0;
        double totalSLOD = 0.0;
        int SLODsSeen = 0;
        double totalFreq = 0.0;
        int freqsSeen = 0;

        for ( RodVCF rod : rods ) {
            List<Genotype> myGenotypes = rod.getGenotypes();
            for ( Genotype g : myGenotypes ) {
                if ( !(g instanceof VCFGenotypeCall) )
                    throw new StingException("Expected VCFGenotypeCall object but instead saw " + g.getClass().getSimpleName());

                // set the name to be the new uniquified name and add it to the list of genotypes
                VCFGenotypeCall call = (VCFGenotypeCall)g;
                call.setSampleName(rodNamesToSampleNames.get(new Pair<String, String>(rod.getName(), call.getSampleName())));
                if ( params.getPosition() < 1 )
                    params.setLocations(call.getLocation(), call.getReference());
                params.addGenotypeRecord(createVCFGenotypeRecord(params, call));
                totalReadDepth += call.getReadCount();
            }

            // set the overall confidence to be the max entry we see
            double confidence = 10.0 * rod.getNegLog10PError();
            if ( confidence > maxConfidence )
                maxConfidence = confidence;

            if ( rod.hasNonRefAlleleFrequency() ) {
                totalFreq += rod.getNonRefAlleleFrequency();
                freqsSeen++;
            }

            if ( rod.hasStrandBias() ) {
                totalSLOD += rod.getStrandBias();
                SLODsSeen++;
            }
        }

        Map<String, String> infoFields = new HashMap<String, String>();
        infoFields.put("DP", String.format("%d", totalReadDepth));
        infoFields.put("NS", String.valueOf(params.getGenotypesRecords().size()));

        // set the overall strand bias and allele frequency to be the average of all entries we've seen
        if ( SLODsSeen > 0 )
            infoFields.put("SB", String.format("%.2f", (totalSLOD/(double)SLODsSeen)));
        if ( freqsSeen > 0 )
            infoFields.put("AF", String.format("%.2f", (totalFreq/(double)freqsSeen)));

        return new VCFRecord(params.getReferenceBase(),
                params.getContig(),
                params.getPosition(),
                ".",
                params.getAlternateBases(),
                maxConfidence,
                "0",
                infoFields,
                params.getFormatString(),
                params.getGenotypesRecords());
    }

    /**
     * create the VCF genotype record
     *
     * @param params the VCF parameters object
     * @param gtype  the genotype
     *
     * @return a VCFGenotypeRecord
     */
    public static VCFGenotypeRecord createVCFGenotypeRecord(VCFParameters params, VCFGenotypeCall gtype) {
        Map<String, String> map = new HashMap<String, String>();

        // calculate the RMS mapping qualities and the read depth
        int readDepth = gtype.getReadCount();
        map.put("RD", String.valueOf(readDepth));
        params.addFormatItem("RD");
        double qual = gtype.getNegLog10PError();
        map.put("GQ", String.format("%.2f", qual));
        params.addFormatItem("GQ");

        List<VCFGenotypeEncoding> alleles = createAlleleArray(gtype);
        for (VCFGenotypeEncoding allele : alleles) {
            params.addAlternateBase(allele);
        }

        VCFGenotypeRecord record = new VCFGenotypeRecord(gtype.getSampleName(),
                                                         alleles,
                                                         VCFGenotypeRecord.PHASE.UNPHASED,
                                                         map);
        return record;
    }

    /**
     * create the allele array?
     *
     * @param gtype the gentoype object
     *
     * @return a list of string representing the string array of alleles
     */
    private static List<VCFGenotypeEncoding> createAlleleArray(Genotype gtype) {
        List<VCFGenotypeEncoding> alleles = new ArrayList<VCFGenotypeEncoding>();
        for (char allele : gtype.getBases().toCharArray()) {
            alleles.add(new VCFGenotypeEncoding(String.valueOf(allele)));
        }
        return alleles;
    }
}