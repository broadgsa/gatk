package org.broadinstitute.sting.utils.genotype.vcf;

import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.OutputStream;
import java.util.*;


/**
 * @author aaron, ebanks
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

    // allowed genotype format strings
    private List<String> allowedGenotypeFormatStrings = null;

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

        // set up the header fields
        Set<VCFHeaderLine> hInfo = new TreeSet<VCFHeaderLine>();
        hInfo.add(new VCFHeaderLine(VCFHeaderVersion.VCF3_3.getFormatString(), VCFHeaderVersion.VCF3_3.getVersionString()));

        // set up the allowed genotype format fields
        if ( headerInfo != null ) {
            for ( VCFHeaderLine field : headerInfo ) {
                hInfo.add(field);
                if ( field instanceof VCFFormatHeaderLine) {
                    if ( allowedGenotypeFormatStrings == null )
                        allowedGenotypeFormatStrings = new ArrayList<String>();
                    allowedGenotypeFormatStrings.add(((VCFFormatHeaderLine)field).getName());
                }
            }
        }

        // set up the sample names
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
     * @param refAllele currently has to be a single character representing the reference base (the base
     * immediately preceding the event in case of indels)
     */
    public void addCall(VariantContext vc, String refAllele) {
        if ( mHeader == null )
            throw new IllegalStateException("The VCF Header must be written before records can be added");

        VCFRecord call = VariantContextAdaptors.toVCF(vc, (byte)refAllele.charAt(0), allowedGenotypeFormatStrings, false, false);

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
