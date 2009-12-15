package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.*;
import java.util.Map.Entry;

/**
 * A set of static utility methods for common operations on VCF files/records.
 */
public class VCFUtils {
    /**
     * Constructor access disallowed...static utility methods only!
     */
    private VCFUtils() { }


    /**
     * Gets the header fields from all VCF rods input by the user
     *
     * @param toolkit    GATK engine
     *
     * @return a set of all fields
     */
    public static Set<VCFHeaderLine> getHeaderFields(GenomeAnalysisEngine toolkit) {

        // keep a map of sample name to occurrences encountered
        TreeSet<VCFHeaderLine> fields = new TreeSet<VCFHeaderLine>();

        // iterate to get all of the sample names
        List<ReferenceOrderedDataSource> dataSources = toolkit.getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            ReferenceOrderedData rod = source.getReferenceOrderedData();
            if ( rod.getType().equals(RodVCF.class) ) {
                VCFReader reader = new VCFReader(rod.getFile());
                fields.addAll(reader.getHeader().getMetaData());
                reader.close();
            }
        }

        return fields;
    }

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

        // keep a map of sample name to occurrences encountered
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

        // how many occurrences have we seen so far?
        Integer occurrences = sampleOverlapMap.get(newSample);

        // if this is the first one, just add it to the list of samples
        if ( occurrences == null ) {
            samples.add(newSample);
            rodNamesToSampleNames.put(new Pair<String, String>(rodName, newSample), newSample);
            sampleOverlapMap.put(newSample, 1);
        }

        // if it's already been seen multiple times, give it a unique suffix and increment the value
        else if ( occurrences >= 2 ) {
            String uniqueName = newSample + "." + rodName;
            samples.add(uniqueName);
            rodNamesToSampleNames.put(new Pair<String, String>(rodName, newSample), uniqueName);
            sampleOverlapMap.put(newSample, occurrences + 1);
        }

        // if this is the second occurrence of the sample name, uniquify both of them
        else { // occurrences == 2

            // remove the 1st occurrence, uniquify it, and add it back
            samples.remove(newSample);
            String uniqueName1 = null;
            for ( Entry<Pair<String, String>, String> entry : rodNamesToSampleNames.entrySet() ) {
                if ( entry.getValue().equals(newSample) ) {
                    uniqueName1 = newSample + "." + entry.getKey().first;
                    entry.setValue(uniqueName1);
                    break;
                }
            }
            samples.add(uniqueName1);

            // add the second one
            String uniqueName2 = newSample + "." + rodName;
            samples.add(uniqueName2);
            rodNamesToSampleNames.put(new Pair<String, String>(rodName, newSample), uniqueName2);

            sampleOverlapMap.put(newSample, 2);
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
        params.addFormatItem(VCFGenotypeRecord.GENOTYPE_KEY);

        // keep track of the locus specific data so we can merge them intelligently
        int totalReadDepth = 0;
        double maxConfidence = 0.0;
        double totalSLOD = 0.0;
        int SLODsSeen = 0;
        double totalFreq = 0.0;
        int freqsSeen = 0;

        for ( RodVCF rod : rods ) {
            List<VCFGenotypeRecord> myGenotypes = rod.getVCFGenotypeRecords();
            for ( VCFGenotypeRecord call : myGenotypes ) {
                // set the name to be the new uniquified name and add it to the list of genotypes
                call.setSampleName(rodNamesToSampleNames.get(new Pair<String, String>(rod.getName(), call.getSampleName())));
                if ( params.getPosition() < 1 )
                    params.setLocations(rod.getLocation(), call.getReference());
                params.addGenotypeRecord(createVCFGenotypeRecord(params, call, rod.mCurrentRecord));
                int depth = call.getReadCount();
                if ( depth > 0 )
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
        infoFields.put(VCFRecord.DEPTH_KEY, String.format("%d", totalReadDepth));
        infoFields.put(VCFRecord.SAMPLE_NUMBER_KEY, String.valueOf(params.getGenotypesRecords().size()));

        // set the overall strand bias and allele frequency to be the average of all entries we've seen
        if ( SLODsSeen > 0 )
            infoFields.put(VCFRecord.STRAND_BIAS_KEY, String.format("%.2f", (totalSLOD/(double)SLODsSeen)));
        if ( freqsSeen > 0 )
            infoFields.put(VCFRecord.ALLELE_FREQUENCY_KEY, String.format("%.2f", (totalFreq/(double)freqsSeen)));

        // TODO -- "." and "0" are wrong -- need to use values from the records

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
     * @param vcfrecord  the VCF record
     *
     * @return a VCFGenotypeRecord
     */
    public static VCFGenotypeRecord createVCFGenotypeRecord(VCFParameters params, VCFGenotypeRecord gtype, VCFRecord vcfrecord) {
        Map<String, String> map = new HashMap<String, String>();

        // calculate the genotype quality and the read depth
        map.put(VCFGenotypeRecord.DEPTH_KEY, String.valueOf(gtype.getReadCount()));
        params.addFormatItem(VCFGenotypeRecord.DEPTH_KEY);
        double qual = Math.min(10.0 * gtype.getNegLog10PError(), VCFGenotypeRecord.MAX_QUAL_VALUE);
        map.put(VCFGenotypeRecord.GENOTYPE_QUALITY_KEY, String.format("%.2f", qual));
        params.addFormatItem(VCFGenotypeRecord.GENOTYPE_QUALITY_KEY);

        List<VCFGenotypeEncoding> alleles = createAlleleArray(gtype);
        for (VCFGenotypeEncoding allele : alleles) {
            params.addAlternateBase(allele);
        }

        VCFGenotypeRecord record = new VCFGenotypeRecord(gtype.getSampleName(),
                                                         alleles,
                                                         VCFGenotypeRecord.PHASE.UNPHASED,
                                                         map);
        record.setVCFRecord(vcfrecord);
        return record;
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
        map.put(VCFGenotypeRecord.DEPTH_KEY, String.valueOf(gtype.getReadCount()));
        params.addFormatItem(VCFGenotypeRecord.DEPTH_KEY);
        double qual = Math.min(10.0 * gtype.getNegLog10PError(), VCFGenotypeRecord.MAX_QUAL_VALUE);
        map.put(VCFGenotypeRecord.GENOTYPE_QUALITY_KEY, String.format("%.2f", qual));
        params.addFormatItem(VCFGenotypeRecord.GENOTYPE_QUALITY_KEY);

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