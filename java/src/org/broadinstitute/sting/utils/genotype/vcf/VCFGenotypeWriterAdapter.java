package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.utils.genotype.*;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.OutputStream;
import java.util.*;


/**
 * @author aaron
 *         <p/>
 *         Class VCFGenotypeWriterAdapter
 *         <p/>
 *         Adapt the VCF writter to the genotype output system
 */
public class VCFGenotypeWriterAdapter implements VCFGenotypeWriter {
    // our VCF objects
    private VCFWriter mWriter = null;
    private VCFHeader mHeader = null;
    private final Set<String> mSampleNames = new LinkedHashSet<String>();

    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(VCFGenotypeWriterAdapter.class);


    public VCFGenotypeWriterAdapter(File writeTo) {
        if (writeTo == null) throw new RuntimeException("VCF output file must not be null");
        mWriter = new VCFWriter(writeTo);
    }

    public VCFGenotypeWriterAdapter(OutputStream writeTo) {
        if (writeTo == null) throw new RuntimeException("VCF output stream must not be null");
        mWriter = new VCFWriter(writeTo);
    }

    /**
     * initialize this VCF header
     *
     * @param sampleNames  the sample names
     * @param headerInfo  the optional header fields
     */
    @Override
    public void writeHeader(Set<String> sampleNames, Set<VCFHeaderLine> headerInfo) {
        mSampleNames.addAll(sampleNames);

        // setup the header fields
        Set<VCFHeaderLine> hInfo = new TreeSet<VCFHeaderLine>();
        hInfo.add(new VCFHeaderLine(VCFHeader.FILE_FORMAT_KEY, VCFHeader.VCF_VERSION));
        hInfo.addAll(headerInfo);
        
        // setup the sample names
        mHeader = new VCFHeader(hInfo, mSampleNames);
        mWriter.writeHeader(mHeader);
    }

    /**
     * Add a genotype, given a genotype locus
     *
     * @param call the locus to add
     */
    public void addGenotypeCall(Genotype call) {
        throw new UnsupportedOperationException("VCF calls require locusdata; use the addMultiSampleCall method instead");
    }

    /**
     * add a no call to the genotype file, if supported.
     *
     * @param position the position to add the no call at
     */
    public void addNoCall(int position) {
        throw new UnsupportedOperationException("We don't currently support no-calls in VCF");
    }

    /** finish writing, closing any open files. */
    public void close() {
        mWriter.close();
    }

    /**
     * add a multi-sample call if we support it
     *
     * @param genotypes the list of genotypes
     */
    public void addMultiSampleCall(List<Genotype> genotypes, VariationCall locusdata) {
        if ( mHeader == null )
            throw new IllegalStateException("The VCF Header must be written before records can be added");

        if ( locusdata != null && !(locusdata instanceof VCFVariationCall) )
            throw new IllegalArgumentException("Only VCFVariationCall objects should be passed in to the VCF writers");

        VCFParameters params = new VCFParameters();
        params.addFormatItem(VCFGenotypeRecord.GENOTYPE_KEY);

        // get the location and reference
        if ( genotypes.size() == 0 ) {
            if ( locusdata == null )
                throw new IllegalArgumentException("Unable to parse out the current location: genotype array must contain at least one entry or have locusdata");

            params.setLocations(locusdata.getLocation(), locusdata.getReference().charAt(0));

            // if there is no genotype data, we'll also need to set an alternate allele
            if ( locusdata.isBiallelic() && locusdata.isSNP() )
                params.addAlternateBase(new VCFGenotypeEncoding(locusdata.getAlternateAlleleList().get(0)));
        } else {
            params.setLocations(genotypes.get(0).getLocation(), genotypes.get(0).getReference());
        }

        Map<String, VCFGenotypeCall> genotypeMap = genotypeListToSampleNameMap(genotypes);

        for (String name : mHeader.getGenotypeSamples()) {
            if (genotypeMap.containsKey(name)) {
                Genotype gtype = genotypeMap.get(name);
                VCFGenotypeRecord record = VCFUtils.createVCFGenotypeRecord(params, (VCFGenotypeCall)gtype);
                params.addGenotypeRecord(record);
                genotypeMap.remove(name);
            } else {
                VCFGenotypeRecord record = createNoCallRecord(name);
                params.addGenotypeRecord(record);
            }
        }

        if (genotypeMap.size() > 0) {
            for (String name : genotypeMap.keySet())
                logger.fatal("Genotype " + name + " was present in the VCFHeader");
            throw new IllegalArgumentException("Genotype array passed to VCFGenotypeWriterAdapter contained Genotypes not in the VCF header");
        }

        // info fields
        Map<String, String> infoFields = getInfoFields((VCFVariationCall)locusdata, params);

        // q-score
        double qual = (locusdata == null) ? 0 : ((VCFVariationCall)locusdata).getConfidence();
        // min Q-score is zero
        qual = Math.max(qual, 0);

        // dbsnp id
        String dbSnpID = null;
        if ( locusdata != null )
            dbSnpID = ((VCFVariationCall)locusdata).getID();

        VCFRecord vcfRecord = new VCFRecord(params.getReferenceBase(),
                                            params.getContig(),
                                            params.getPosition(),
                                            (dbSnpID == null ? VCFRecord.EMPTY_ID_FIELD : dbSnpID),
                                            params.getAlternateBases(),
                                            qual,
                                            VCFRecord.UNFILTERED,
                                            infoFields,
                                            params.getFormatString(),
                                            params.getGenotypesRecords());

        mWriter.addRecord(vcfRecord);
    }

    public void addRecord(VCFRecord vcfRecord) {
        mWriter.addRecord(vcfRecord);
    }

    /**
     * get the information fields of the VCF record, given the meta data and parameters
     *
     * @param locusdata the metadata associated with this multi sample call
     * @param params   the parameters
     *
     * @return a mapping of info field to value
     */
    private static Map<String, String> getInfoFields(VCFVariationCall locusdata, VCFParameters params) {
        Map<String, String> infoFields = new HashMap<String, String>();
        if ( locusdata != null ) {
            if ( locusdata.getSLOD() != null )
                infoFields.put(VCFRecord.STRAND_BIAS_KEY, String.format("%.2f", locusdata.getSLOD()));
            if ( locusdata.hasNonRefAlleleFrequency() )
                infoFields.put(VCFRecord.ALLELE_FREQUENCY_KEY, String.format("%.2f", locusdata.getNonRefAlleleFrequency()));
            Map<String, String> otherFields = locusdata.getFields();
            if ( otherFields != null ) {
                infoFields.putAll(otherFields);
            }
        }
        infoFields.put(VCFRecord.SAMPLE_NUMBER_KEY, String.valueOf(params.getGenotypesRecords().size()));
        return infoFields;
    }

    /**
     * create a no call record
     *
     * @param sampleName the sample name
     *
     * @return a VCFGenotypeRecord for the no call situation
     */
    private VCFGenotypeRecord createNoCallRecord(String sampleName) {
        Map<String, String> map = new HashMap<String, String>();


        List<VCFGenotypeEncoding> alleles = new ArrayList<VCFGenotypeEncoding>();
        alleles.add(new VCFGenotypeEncoding(VCFGenotypeRecord.EMPTY_ALLELE));
        alleles.add(new VCFGenotypeEncoding(VCFGenotypeRecord.EMPTY_ALLELE));

        VCFGenotypeRecord record = new VCFGenotypeRecord(sampleName,
                                                         alleles,
                                                         VCFGenotypeRecord.PHASE.UNPHASED,
                                                         map);
        return record;
    }

    /** @return true if we support multisample, false otherwise */
    public boolean supportsMultiSample() {
        return true;
    }

    /**
     * create a genotype mapping from a list and their sample names
     * while we're at it, checks that all genotypes are VCF-based
     *
     * @param list a list of genotype samples
     *
     * @return a mapping of the sample name to genotype fields
     */
    private static Map<String, VCFGenotypeCall> genotypeListToSampleNameMap(List<Genotype> list) {
        Map<String, VCFGenotypeCall> map = new HashMap<String, VCFGenotypeCall>();
        for (Genotype rec : list) {
            if ( !(rec instanceof VCFGenotypeCall) )
                throw new IllegalArgumentException("Only VCFGenotypeCalls should be passed in to the VCF writers");
            map.put(((VCFGenotypeCall) rec).getSampleName(), (VCFGenotypeCall) rec);
        }
        return map;
    }

}
