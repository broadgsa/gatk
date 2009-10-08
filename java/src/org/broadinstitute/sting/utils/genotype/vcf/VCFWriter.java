package org.broadinstitute.sting.utils.genotype.vcf;


import java.io.*;
import java.nio.charset.Charset;

/** this class writers VCF files */
public class VCFWriter {

    // the VCF header we're storing
    private VCFHeader mHeader;

    // the print stream we're writting to
    BufferedWriter mWriter;
    private final String FIELD_SEPERATOR = "\t";

    /**
     * create a VCF writer, given a VCF header and a file to write to
     *
     * @param header   the VCF header
     * @param location the file location to write to
     */
    public VCFWriter(VCFHeader header, File location) {
        this.mHeader = header;
        Charset utf8 = Charset.forName("UTF-8");
        try {
            mWriter = new BufferedWriter(
                    new OutputStreamWriter(
                            new FileOutputStream(location),
                            utf8));
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Unable to create VCF file: " + location, e);
        }
        try {
            // write the header meta-data out
            for (String metadata : header.getMetaData().keySet()) {
                mWriter.write(VCFHeader.METADATA_INDICATOR + metadata + "=" + header.getMetaData().get(metadata) + "\n");
            }
            // write out the column line
            StringBuilder b = new StringBuilder();
            b.append(VCFHeader.HEADER_INDICATOR);
            for (VCFHeader.HEADER_FIELDS field : header.getHeaderFields()) b.append(field + FIELD_SEPERATOR);
            if (header.hasGenotypingData()) {
                b.append("FORMAT" + FIELD_SEPERATOR);
                for (String field : header.getGenotypeSamples()) b.append(field + FIELD_SEPERATOR);
                mWriter.write(b.toString() + "\n");
            }
        }
        catch (IOException e) {
            throw new RuntimeException("IOException writing the VCF header", e);
        }
    }

    /**
     * output a record to the VCF file
     *
     * @param record the record to output
     */
    public void addRecord(VCFRecord record) {
        if (record.getColumnCount() != mHeader.getGenotypeSamples().size() + mHeader.getHeaderFields().size()) {
            throw new RuntimeException("Record has " + record.getColumnCount() +
                    " columns, when is should have " + mHeader.getColumnCount());
        }
        String vcfString = record.toString();
        try {
            mWriter.write(vcfString + "\n");
        } catch (IOException e) {
            throw new RuntimeException("Unable to write the VCF object to a file");
        }

    }



    /** attempt to close the VCF file */
    public void close() {
        try {
            mWriter.flush();
            mWriter.close();
        } catch (IOException e) {
            throw new RuntimeException("Unable to close VCFFile");
        }
    }

}
