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
public class VCFGenotypeWriterAdapter implements GenotypeWriter {
    // our VCF objects
    private VCFWriter mWriter = null;
    private VCFHeader mHeader = null;
    private String mSource;
    private String mReferenceName;
    private boolean mInitialized = false;
    private final Set<String> mSampleNames = new HashSet<String>();
    private final File mFile;
    private final OutputStream mStream;

    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(VCFGenotypeWriterAdapter.class);


    public VCFGenotypeWriterAdapter(String source, String referenceName, File writeTo, Set<String> sampleNames) {
        mReferenceName = referenceName;
        mSource = source;
        mFile = writeTo;
        if (mFile == null) throw new RuntimeException("VCF output file must not be null");
        mStream = null;
        mSampleNames.addAll(sampleNames);
    }

    public VCFGenotypeWriterAdapter(String source, String referenceName, OutputStream writeTo, Set<String> sampleNames) {
        mReferenceName = referenceName;
        mSource = source;
        mFile = null;
        mStream = writeTo;
        if (mStream == null) throw new RuntimeException("VCF output stream must not be null");
        mSampleNames.addAll(sampleNames);

    }

    /**
     * initialize this VCF writer
     *
     * @param genotypes the genotypes
     * @param file      the file location to write to
     */
    private void lazyInitialize(List<Genotype> genotypes, File file, OutputStream stream) {
        Map<String, String> hInfo = new HashMap<String, String>();

        // setup the header fields
        hInfo.put("format", "VCRv3.2");
        hInfo.put("source", mSource);
        hInfo.put("reference", mReferenceName);

        // setup the sample names
        mHeader = new VCFHeader(hInfo, mSampleNames);
        if (mFile == null)
            mWriter = new VCFWriter(mHeader, stream);
        else
            mWriter = new VCFWriter(mHeader, file);
        mInitialized = true;
    }

    /**
     * get the samples names from genotype objects
     *
     * @param genotypes the genotype list
     *
     * @return a list of strings representing the sample names
     */
    private static List<String> getSampleNames(List<Genotype> genotypes) {
        List<String> strings = new ArrayList<String>();
        for (Genotype genotype : genotypes) {
            if (!(genotype instanceof SampleBacked))
                throw new IllegalArgumentException("Genotypes passed to VCF must be backed by SampledBacked interface");
            strings.add(((SampleBacked) genotype).getSampleName());
        }
        return strings;
    }

    /**
     * Add a genotype, given a genotype locus
     *
     * @param call the locus to add
     */
    @Override
    public void addGenotypeCall(Genotype call) {
        addMultiSampleCall(Arrays.asList(call), null);
    }

    /**
     * add a no call to the genotype file, if supported.
     *
     * @param position the position to add the no call at
     */
    @Override
    public void addNoCall(int position) {
        throw new UnsupportedOperationException("We don't currently support no-calls in VCF");
    }

    /** finish writing, closing any open files. */
    @Override
    public void close() {
        if (mInitialized)
            mWriter.close();
    }

    /**
     * add a multi-sample call if we support it
     *
     * @param genotypes the list of genotypes, that are backed by sample information
     */
    @Override
    public void addMultiSampleCall(List<Genotype> genotypes, GenotypeMetaData metadata) {
        if (!mInitialized)
            lazyInitialize(genotypes, mFile, mStream);


        VCFParameters params = new VCFParameters();
        params.addFormatItem("GT");

        // mapping of our sample names to genotypes
        if (genotypes.size() < 1) {
            throw new IllegalArgumentException("Unable to parse out the current location: genotype array must at least contain one entry");
        }

        // get the location and reference
        params.setLocations(genotypes.get(0).getLocation(), genotypes.get(0).getReference());

        Map<String, Genotype> genotypeMap = genotypeListToSampleNameMap(genotypes);

        for (String name : mHeader.getGenotypeSamples()) {
            if (genotypeMap.containsKey(name)) {
                Genotype gtype = genotypeMap.get(name);
                VCFGenotypeRecord record = createVCFGenotypeRecord(params, gtype);
                params.addGenotypeRecord(record);
                genotypeMap.remove(name);
            } else {
                VCFGenotypeRecord record = createNoCallRecord(params, name);
                params.addGenotypeRecord(record);
            }
        }

        if (genotypeMap.size() > 0) {
            for (String name : genotypeMap.keySet())
                logger.fatal("Genotype " + name + " was present in the VCFHeader");
            throw new IllegalArgumentException("Genotype array passed to VCFGenotypeWriterAdapter contained Genotypes not in the VCF header");
        }

        Map<String, String> infoFields = getInfoFields(metadata, params);

        double qual = (metadata == null) ? 0 : (metadata.getLOD()) * 10;

        /**
         * TODO: Eric fix the next line when our LOD scores are 0->Inf based instead
         * of -3 to Inf based.
         */
        if (qual < 0.0) {
            qual = 0.0;
        }

        VCFRecord vcfRecord = new VCFRecord(params.getReferenceBase(),
                                            params.getContig(),
                                            params.getPosition(),
                                            ".",
                                            params.getAlternateBases(),
                                            qual,
                                            ".",
                                            infoFields,
                                            params.getFormatString(),
                                            params.getGenotypesRecords());

        mWriter.addRecord(vcfRecord);
    }

    /**
     * get the information fields of the VCF record, given the meta data and parameters
     *
     * @param metadata the metadata associated with this multi sample call
     * @param params   the parameters
     *
     * @return a mapping of info field to value
     */
    private Map<String, String> getInfoFields(GenotypeMetaData metadata, VCFParameters params) {
        Map<String, String> infoFields = new HashMap<String, String>();
        if (metadata != null) {
            infoFields.put("SB", String.format("%.2f", metadata.getSLOD()));
            infoFields.put("AF", String.format("%.2f", metadata.getAlleleFrequency()));
        }
        infoFields.put("NS", String.valueOf(params.getGenotypesRecords().size()));
        return infoFields;
    }

    /**
     * create the VCF genotype record
     *
     * @param params the VCF parameters object
     * @param gtype  the genotype
     *
     * @return a VCFGenotypeRecord
     */
    private VCFGenotypeRecord createVCFGenotypeRecord(VCFParameters params, Genotype gtype) {
        Map<String, String> map = new HashMap<String, String>();
        if (!(gtype instanceof SampleBacked)) {
            throw new IllegalArgumentException("Genotypes passed to VCF must be backed by SampledBacked interface");
        }

        // calculate the RMS mapping qualities and the read depth
        if (gtype instanceof ReadBacked) {
            int readDepth = ((ReadBacked) gtype).getReadCount();
            map.put("RD", String.valueOf(readDepth));
            params.addFormatItem("RD");
        }
        double qual = gtype.getNegLog10PError();
        map.put("GQ", String.format("%.2f", qual));
        params.addFormatItem("GQ");

        List<VCFGenotypeEncoding> alleles = createAlleleArray(gtype);
        for (VCFGenotypeEncoding allele : alleles) {
            params.addAlternateBase(allele);
        }

        VCFGenotypeRecord record = new VCFGenotypeRecord(((SampleBacked) gtype).getSampleName(),
                                                         alleles,
                                                         VCFGenotypeRecord.PHASE.UNPHASED,
                                                         map);
        return record;
    }

    /**
     * create a no call record
     *
     * @param params     the VCF parameters object
     * @param sampleName the sample name
     *
     * @return a VCFGenotypeRecord for the no call situation
     */
    private VCFGenotypeRecord createNoCallRecord(VCFParameters params, String sampleName) {
        Map<String, String> map = new HashMap<String, String>();


        List<VCFGenotypeEncoding> alleles = new ArrayList<VCFGenotypeEncoding>();
        alleles.add(new VCFGenotypeEncoding(VCFGenotypeRecord.EMPTY_GENOTYPE));
        alleles.add(new VCFGenotypeEncoding(VCFGenotypeRecord.EMPTY_GENOTYPE));

        VCFGenotypeRecord record = new VCFGenotypeRecord(sampleName,
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
    private List<VCFGenotypeEncoding> createAlleleArray(Genotype gtype) {
        List<VCFGenotypeEncoding> alleles = new ArrayList<VCFGenotypeEncoding>();
        for (char allele : gtype.getBases().toCharArray()) {
            alleles.add(new VCFGenotypeEncoding(String.valueOf(allele)));
        }
        return alleles;
    }

    /** @return true if we support multisample, false otherwise */
    @Override
    public boolean supportsMultiSample() {
        return true;
    }

    /**
     * create a genotype mapping from a list and their sample names
     *
     * @param list a list of genotype samples
     *
     * @return a mapping of the sample name to genotype fields
     */
    private static Map<String, Genotype> genotypeListToSampleNameMap(List<Genotype> list) {
        Map<String, Genotype> map = new HashMap<String, Genotype>();
        for (Genotype rec : list) {
            if (!(rec instanceof SampleBacked))
                throw new IllegalArgumentException("Genotype must be backed by sample information");
            map.put(((SampleBacked) rec).getSampleName(), rec);
        }
        return map;
    }

}
