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
     * @param file      the file location to write to
     * @param stream    the output stream
     */
    private void lazyInitialize(File file, OutputStream stream) {
        Map<String, String> hInfo = new HashMap<String, String>();

        // setup the header fields
        hInfo.put("format", VCFWriter.VERSION);
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
            if (!(genotype instanceof VCFGenotypeCall))
                throw new IllegalArgumentException("Genotypes passed to VCF must be backed by SampledBacked interface");
            strings.add(((VCFGenotypeCall) genotype).getSampleName());
        }
        return strings;
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
        if (mInitialized)
            mWriter.close();
    }

    /**
     * add a multi-sample call if we support it
     *
     * @param genotypes the list of genotypes
     */
    public void addMultiSampleCall(List<Genotype> genotypes, GenotypeLocusData locusdata) {
        if (!mInitialized)
            lazyInitialize(mFile, mStream);

        if ( locusdata != null && !(locusdata instanceof VCFGenotypeLocusData) )
            throw new IllegalArgumentException("Only VCFGenotypeLocusData objects should be passed in to the VCF writers");

        VCFParameters params = new VCFParameters();
        params.addFormatItem("GT");

        // mapping of our sample names to genotypes
        if (genotypes.size() < 1) {
            throw new IllegalArgumentException("Unable to parse out the current location: genotype array must contain at least one entry");
        }

        // get the location and reference
        params.setLocations(genotypes.get(0).getLocation(), genotypes.get(0).getReference());

        Map<String, VCFGenotypeCall> genotypeMap = genotypeListToSampleNameMap(genotypes);

        // keep track of the total read depth
        int totalReadDepth = 0;

        for (String name : mHeader.getGenotypeSamples()) {
            if (genotypeMap.containsKey(name)) {
                Genotype gtype = genotypeMap.get(name);
                VCFGenotypeRecord record = VCFUtils.createVCFGenotypeRecord(params, (VCFGenotypeCall)gtype);
                totalReadDepth += ((VCFGenotypeCall)gtype).getReadCount();
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

        Map<String, String> infoFields = getInfoFields((VCFGenotypeLocusData)locusdata, params);
        infoFields.put("DP", String.format("%d", totalReadDepth));

        double qual = (locusdata == null) ? 0 : ((VCFGenotypeLocusData)locusdata).getConfidence();
        // maintain 0-99 based Q-scores
        qual = Math.min(qual, 99);
        qual = Math.max(qual, 0);

        VCFRecord vcfRecord = new VCFRecord(params.getReferenceBase(),
                                            params.getContig(),
                                            params.getPosition(),
                                            ".",
                                            params.getAlternateBases(),
                                            qual,
                                            "0",
                                            infoFields,
                                            params.getFormatString(),
                                            params.getGenotypesRecords());

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
    private static Map<String, String> getInfoFields(VCFGenotypeLocusData locusdata, VCFParameters params) {
        Map<String, String> infoFields = new HashMap<String, String>();
        if ( locusdata != null ) {
            infoFields.put("SB", String.format("%.2f", locusdata.getSLOD()));
            infoFields.put("AF", String.format("%.2f", locusdata.getAlleleFrequency()));
        }
        infoFields.put("NS", String.valueOf(params.getGenotypesRecords().size()));
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
        alleles.add(new VCFGenotypeEncoding(VCFGenotypeRecord.EMPTY_GENOTYPE));
        alleles.add(new VCFGenotypeEncoding(VCFGenotypeRecord.EMPTY_GENOTYPE));

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
