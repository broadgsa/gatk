package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.Utils;

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

        // keep track of the data so we can merge them intelligently
        double maxConfidence = 0.0;
        String id = null;
        Map<String, String> infoFields = new HashMap<String, String>();
        List<String> filters = new ArrayList<String>();

        for ( RodVCF rod : rods ) {
            List<VCFGenotypeRecord> myGenotypes = rod.getVCFGenotypeRecords();
            for ( VCFGenotypeRecord call : myGenotypes ) {
                // set the name to be the new uniquified name and add it to the list of genotypes
                call.setSampleName(rodNamesToSampleNames.get(new Pair<String, String>(rod.getName(), call.getSampleName())));
                if ( params.getPosition() < 1 )
                    params.setLocations(rod.getLocation(), call.getReference());
                params.addGenotypeRecord(createVCFGenotypeRecord(params, call, rod.mCurrentRecord));
            }

            // set the overall confidence to be the max entry we see
            double confidence = 10.0 * rod.getNegLog10PError();
            if ( confidence > maxConfidence )
                maxConfidence = confidence;

            if ( rod.getID() != null )
                id = rod.getID();

            if ( rod.isFiltered() )
                filters.add(rod.getFilterString());

            // just take the last value we see for a given key
            infoFields.putAll(rod.getInfoValues());
        }

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
        for ( Map.Entry<String, String> entry : gtype.getFields().entrySet() ) {
            record.setField(entry.getKey(), entry.getValue());
            params.addFormatItem(entry.getKey());
        }

        record.setVCFRecord(vcfrecord);
        return record;
    }

    /**
     * create the allele array?
     *
     * @param gtype the gentoype object
     *
     * @return a list of string representing the string array of alleles
     */
    private static List<VCFGenotypeEncoding> createAlleleArray(VCFGenotypeRecord gtype) {
        List<VCFGenotypeEncoding> alleles = new ArrayList<VCFGenotypeEncoding>();
        for (char allele : gtype.getBases().toCharArray()) {
            alleles.add(new VCFGenotypeEncoding(String.valueOf(allele)));
        }
        return alleles;
    }
}