package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.*;

/**
 * A set of static utility methods for common operations on VCF files/records.
 */
public class VCFUtils {
    /**
     * Constructor access disallowed...static utility methods only!
     */
    private VCFUtils() { }

    public static Set<ReferenceOrderedData> getRodVCFs(GenomeAnalysisEngine toolkit) {
        Set<ReferenceOrderedData> vcfs = new HashSet<ReferenceOrderedData>();

        for ( ReferenceOrderedDataSource source : toolkit.getRodDataSources() ) {
            ReferenceOrderedData rod = source.getReferenceOrderedData();
            if ( rod.getType().equals(RodVCF.class) ) {
                vcfs.add(rod);
            }
        }

        return vcfs;
    }

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
        String id = null;
        List<String> filters = new ArrayList<String>();

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

            if ( !rod.isReference() && rod.hasNonRefAlleleFrequency() ) {
                totalFreq += rod.getNonRefAlleleFrequency();
                freqsSeen++;
            }

            if ( rod.hasStrandBias() ) {
                totalSLOD += rod.getStrandBias();
                SLODsSeen++;
            }

            if ( rod.getID() != null )
                id = rod.getID();

            if ( rod.isFiltered() )
                filters.add(rod.getFilterString());
        }

        Map<String, String> infoFields = new HashMap<String, String>();
        infoFields.put(VCFRecord.DEPTH_KEY, String.format("%d", totalReadDepth));

        // set the overall strand bias and allele frequency to be the average of all entries we've seen
        if ( SLODsSeen > 0 )
            infoFields.put(VCFRecord.STRAND_BIAS_KEY, String.format("%.2f", (totalSLOD/(double)SLODsSeen)));
        if ( freqsSeen > 0 )
            infoFields.put(VCFRecord.ALLELE_FREQUENCY_KEY, String.format("%.2f", (totalFreq/(double)freqsSeen)));
                
        return new VCFRecord(params.getReferenceBases(),
                params.getContig(),
                params.getPosition(),
                (id != null ? id : "."),
                params.getAlternateBases(),
                maxConfidence,
                filters.size() == 0 ? "0" : Utils.join(";", filters),
                infoFields,
                params.getFormatString(),
                params.getGenotypeRecords());
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

        List<VCFGenotypeEncoding> alleles = createAlleleArray(gtype);
        for (VCFGenotypeEncoding allele : alleles) {
            params.addAlternateBase(allele);
        }

        VCFGenotypeRecord record = new VCFGenotypeRecord(gtype.getSampleName(), alleles, VCFGenotypeRecord.PHASE.UNPHASED);

        // calculate the genotype quality and the read depth
        record.setField(VCFGenotypeRecord.DEPTH_KEY, String.valueOf(gtype.getReadCount()));
        params.addFormatItem(VCFGenotypeRecord.DEPTH_KEY);
        double qual = Math.min(10.0 * gtype.getNegLog10PError(), VCFGenotypeRecord.MAX_QUAL_VALUE);
        record.setField(VCFGenotypeRecord.GENOTYPE_QUALITY_KEY, String.format("%.2f", qual));
        params.addFormatItem(VCFGenotypeRecord.GENOTYPE_QUALITY_KEY);

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

        List<VCFGenotypeEncoding> alleles = createAlleleArray(gtype);
        for (VCFGenotypeEncoding allele : alleles) {
            params.addAlternateBase(allele);
        }

        VCFGenotypeRecord record = new VCFGenotypeRecord(gtype.getSampleName(), alleles, VCFGenotypeRecord.PHASE.UNPHASED);

        // calculate the RMS mapping qualities and the read depth
        record.setField(VCFGenotypeRecord.DEPTH_KEY, String.valueOf(gtype.getReadCount()));
        params.addFormatItem(VCFGenotypeRecord.DEPTH_KEY);
        double qual = Math.min(10.0 * gtype.getNegLog10PError(), VCFGenotypeRecord.MAX_QUAL_VALUE);
        record.setField(VCFGenotypeRecord.GENOTYPE_QUALITY_KEY, String.format("%.2f", qual));
        params.addFormatItem(VCFGenotypeRecord.GENOTYPE_QUALITY_KEY);

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