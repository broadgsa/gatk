package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.utils.genotype.*;

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
    private final File mFile;
    private final OutputStream mStream;

    public VCFGenotypeWriterAdapter(String source, String referenceName, File writeTo) {
        mReferenceName = referenceName;
        mSource = source;
        mFile = writeTo;
        if (mFile == null) throw new RuntimeException("VCF output file must not be null");
        mStream = null;
    }

    public VCFGenotypeWriterAdapter(String source, String referenceName, OutputStream writeTo) {
        mReferenceName = referenceName;
        mSource = source;
        mFile = null;
        mStream = writeTo;
        if (mStream == null) throw new RuntimeException("VCF output stream must not be null");

    }

    /**
     * initialize this VCF writer
     *
     * @param genotypes the genotypes
     * @param file      the file location to write to
     */
    private void lazyInitialize(List<Genotype> genotypes, File file, OutputStream stream) {
        Map<String, String> hInfo = new HashMap<String, String>();
        List<String> sampleNames = getSampleNames(genotypes);

        // setup the header fields
        hInfo.put("format", "VCRv3.2");
        hInfo.put("source", mSource);
        hInfo.put("reference", mReferenceName);

        // setup the sample names
        mHeader = new VCFHeader(hInfo, sampleNames);
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

        for (Genotype gtype : genotypes) {
            // setup the parameters
            params.setLocations(gtype.getLocation(), gtype.getReference());

            VCFGenotypeRecord record = createVCFGenotypeRecord(params, gtype);
            params.addGenotypeRecord(record);
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

}
