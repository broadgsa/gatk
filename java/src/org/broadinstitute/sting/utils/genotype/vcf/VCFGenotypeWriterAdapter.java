package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;
import org.broadinstitute.sting.utils.genotype.ReadBacked;
import org.broadinstitute.sting.utils.genotype.SampleBacked;

import java.io.File;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


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
    private final Map<String, String> mSampleNames = new HashMap<String, String>();
    private boolean mInitialized = false;
    private final File mFile;

    public VCFGenotypeWriterAdapter(String source, String referenceName, File writeTo) {
        mReferenceName = referenceName;
        mSource = source;
        mFile = writeTo;
    }


    /**
     * initialize this VCF writer
     *
     * @param genotypes the genotypes
     * @param file      the file location to write to
     */
    private void lazyInitialize(List<Genotype> genotypes, File file) {
        Map<String, String> hInfo = new HashMap<String, String>();
        List<String> sampleNames = getSampleNames(genotypes);

        // setup the header fields
        hInfo.put("format", "VCRv3.2");
        hInfo.put("source", mSource);

        // setup the sample names
        mHeader = new VCFHeader(hInfo, sampleNames);
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
        addMultiSampleCall(Arrays.asList(call));
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
        mWriter.close();
    }

    /**
     * add a multi-sample call if we support it
     *
     * @param genotypes the list of genotypes, that are backed by sample information
     */
    @Override
    public void addMultiSampleCall(List<Genotype> genotypes) {
        if (!mInitialized)
            lazyInitialize(genotypes, mFile);


        VCFParamters params = new VCFParamters();
        params.addFormatItem("GT");

        for (Genotype gtype : genotypes) {
            // setup the parameters
            params.setLocations(gtype.getLocation(), gtype.getReference());

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

            List<String> alleles = new ArrayList<String>();
            for (char allele : gtype.getBases().toCharArray()) {
                alleles.add(String.valueOf(allele));
                params.addAlternateBase(allele);
            }

            VCFGenotypeRecord record = new VCFGenotypeRecord(((SampleBacked) gtype).getSampleName(),
                                                             alleles,
                                                             VCFGenotypeRecord.PHASE.UNPHASED,
                                                             map);
            params.addGenotypeRecord(record);
        }

        Map<String, String> infoFields = new HashMap<String, String>();

        VCFRecord vcfRecord = new VCFRecord(params.getReferenceBase(),
                                            params.getContig(),
                                            params.getPosition(),
                                            ".",
                                            params.getAlternateBases(),
                                            0, /* BETTER VALUE HERE */
                                            ".",
                                            infoFields,
                                            params.getFormatString(),
                                            params.getGenotypesRecords());

        mWriter.addRecord(vcfRecord);
    }

    /** @return true if we support multisample, false otherwise */
    @Override
    public boolean supportsMulitSample() {
        return true;
    }


    /**
     * a helper class, which performs a lot of the safety checks on the parameters
     * we feed to the VCF (like ensuring the same position for each genotype in a call).
     */
    class VCFParamters {
        private char referenceBase = '0';
        private int position = 0;
        private String contig = null;
        private boolean initialized = false;
        private List<VCFGenotypeRecord> genotypesRecord = new ArrayList<VCFGenotypeRecord>();
        private List<String> formatList = new ArrayList<String>();
        private List<String> alternateBases = new ArrayList<String>();

        public void setLocations(GenomeLoc location, char refBase) {
            // if we haven't set it up, we initialize the object
            if (!initialized) {
                initialized = true;
                this.contig = location.getContig();
                this.position = (int)location.getStart();
                if (location.getStart() != location.getStop()) {
                    throw new IllegalArgumentException("The start and stop locations must be the same");
                }
                this.referenceBase = refBase;
            } else {
                if (contig.equals(this.contig))
                    throw new IllegalArgumentException("The contig name has to be the same at a single locus");
                if (position != this.position)
                    throw new IllegalArgumentException("The position has to be the same at a single locus");
                if (refBase != this.referenceBase)
                    throw new IllegalArgumentException("The reference base name has to be the same at a single locus");
            }
        }

        /** @return get the position */
        public int getPosition() {
            return position;
        }

        /** @return get the contig name */
        public String getContig() {
            return contig;
        }

        /** @return get the reference base */
        public char getReferenceBase() {
            return referenceBase;
        }

        public void addGenotypeRecord(VCFGenotypeRecord record) {
            this.genotypesRecord.add(record);
        }

        public void addFormatItem(String item) {
            if (!formatList.contains(item))
                formatList.add(item);
        }

        public void addAlternateBase(char base) {
            if (!alternateBases.contains(String.valueOf(base)) && base != this.getReferenceBase())
                alternateBases.add(String.valueOf(base));
        }

        public List<String> getAlternateBases() {
            return alternateBases;
        }

        public String getFormatString() {
            return Utils.join(";", formatList);
        }

        public List<VCFGenotypeRecord> getGenotypesRecords() {
            return genotypesRecord;
        }
    }
}
