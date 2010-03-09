package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
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

    // our log, which we want to capture anything from this class
    protected static Logger logger = Logger.getLogger(VCFGenotypeWriterAdapter.class);

    // validation stringency
    private VALIDATION_STRINGENCY validationStringency = VALIDATION_STRINGENCY.STRICT;

    // standard genotype format strings
    private static String[] standardGenotypeFormatStrings = { VCFGenotypeRecord.GENOTYPE_KEY,
                                                              VCFGenotypeRecord.DEPTH_KEY,
                                                              VCFGenotypeRecord.GENOTYPE_QUALITY_KEY };

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

    /** finish writing, closing any open files. */
    public void close() {
        mWriter.close();
    }

    /**
     * Add a genotype, given a variant context
     *
     * @param vc  the variant context representing the call to add
     */
    public void addCall(VariantContext vc) {
        if ( mHeader == null )
            throw new IllegalStateException("The VCF Header must be written before records can be added");

        List<String> formatStrings;
        if ( vc.getChromosomeCount() > 0 )
            formatStrings = Arrays.asList(standardGenotypeFormatStrings);
        else
            formatStrings = new ArrayList<String>();
        VCFRecord call = VariantContextAdaptors.toVCF(vc, vc.getReference().toString().charAt(0), formatStrings, false);

        Set<Allele> altAlleles = vc.getAlternateAlleles();
        StringBuffer altAlleleCountString = new StringBuffer();
        for ( Allele allele : altAlleles ) {
            if ( altAlleleCountString.length() > 0 )
                altAlleleCountString.append(",");
            altAlleleCountString.append(vc.getChromosomeCount(allele));
        }
        if ( vc.getChromosomeCount() > 0 ) {
            call.addInfoField(VCFRecord.ALLELE_NUMBER_KEY, String.format("%d", vc.getChromosomeCount()));
            if ( altAlleleCountString.length() > 0 )
                call.addInfoField(VCFRecord.ALLELE_COUNT_KEY, altAlleleCountString.toString());
        }

        mWriter.addRecord(call, validationStringency);
    }

    public void addRecord(VCFRecord vcfRecord) {
        mWriter.addRecord(vcfRecord, validationStringency);
    }

    /**
     * set the validation stringency
     * @param value   validation stringency value
     */
    public void setValidationStringency(VALIDATION_STRINGENCY value) {
        validationStringency = value;
    }
}
