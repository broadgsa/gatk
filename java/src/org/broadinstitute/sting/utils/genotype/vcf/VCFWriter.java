package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.utils.StingException;

import java.io.*;
import java.nio.charset.Charset;

/** this class writers VCF files */
public class VCFWriter {

    // the VCF header we're storing
    private VCFHeader mHeader;

    // the print stream we're writting to
    BufferedWriter mWriter;

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
            throw new StingException("Unable to create VCF file: " + location, e);
        }
        try {

            // write the header meta-data out
            for (String metadata : header.getMetaData().keySet()) {
                mWriter.write(VCFHeader.METADATA_INDICATOR + metadata + "=" + header.getMetaData().get(metadata) + "\n");
            }
            // write out the column line
            StringBuilder b = new StringBuilder();
            b.append(VCFHeader.HEADER_INDICATOR);
            for (VCFHeader.HEADER_FIELDS field : header.getHeaderFields()) b.append(field + "\t");
            for (String field : header.getAuxillaryTags()) b.append(field + "\t");
            mWriter.write(b.toString() + "\n");
        }
        catch (IOException e) {
            throw new StingException("IOException writing the VCF header", e);
        }
    }

    /**
     * output a record to the VCF file
     * @param record the record to output
     */
    public void addRecord(VCFRecord record) {
        if (record.getColumnCount() != mHeader.getAuxillaryTags().size() + mHeader.getHeaderFields().size()) {
            throw new StingException("Record has " + record.getColumnCount() +
                    " columns, when is should have " + (mHeader.getAuxillaryTags().size() +
                    mHeader.getHeaderFields().size()));
        }
        StringBuilder builder = new StringBuilder();
        // first output the required fields in order
        boolean first = true;
        for (VCFHeader.HEADER_FIELDS field : mHeader.getHeaderFields()) {
            if (first) { first = false; builder.append(record.getValue(field)); }
            else builder.append("\t" + record.getValue(field));
        }
        for (String auxTag : mHeader.getAuxillaryTags()) {
            builder.append("\t" + record.getValue(auxTag));
        }
        try {
            mWriter.write(builder.toString() + "\n");
        } catch (IOException e) {
            throw new StingException("Unable to write the VCF object to a file");
        }
    }

    /**
     * attempt to close the VCF file
     */
    public void close() {
        try {
            mWriter.close();
        } catch (IOException e) {
            throw new StingException("Unable to close VCFFile");
        }
    }


}
