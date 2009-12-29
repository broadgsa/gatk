package org.broadinstitute.sting.utils.genotype.vcf;


import java.io.*;
import java.util.TreeSet;

/**
 * this class writers VCF files
 */
public class VCFWriter {

    
    // the VCF header we're storing
    private VCFHeader mHeader = null;

    // the print stream we're writting to
    BufferedWriter mWriter;
    private final String FIELD_SEPERATOR = "\t";

    /**
     * create a VCF writer, given a file to write to
     *
     * @param location the file location to write to
     */
    public VCFWriter(File location) {
        FileOutputStream output;
        try {
            output = new FileOutputStream(location);
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Unable to create VCF file at location: " + location);
        }
        mWriter = new BufferedWriter(new OutputStreamWriter(output));
    }


    /**
     * create a VCF writer, given a stream to write to
     *
     * @param output   the file location to write to
     */
    public VCFWriter(OutputStream output) {
        mWriter = new BufferedWriter(new OutputStreamWriter(output));
    }

    public void writeHeader(VCFHeader header) {
        this.mHeader = header;
        try {
            // the fileformat field needs to be written first
            TreeSet<VCFHeaderLine> nonFormatMetaData = new TreeSet<VCFHeaderLine>();
            for ( VCFHeaderLine line : header.getMetaData() ) {
                if ( line.getKey().equals(VCFHeader.FILE_FORMAT_KEY) ) {
                    mWriter.write(VCFHeader.METADATA_INDICATOR + line.toString() + "\n");
                }
                else if ( line.getKey().equals(VCFHeader.OLD_FILE_FORMAT_KEY) ) {
                    mWriter.write(VCFHeader.METADATA_INDICATOR + VCFHeader.FILE_FORMAT_KEY + line.toString().substring(VCFHeader.OLD_FILE_FORMAT_KEY.length()) + "\n");
                } else {
                    nonFormatMetaData.add(line);
                }
            }

            // write the rest of the header meta-data out
            for ( VCFHeaderLine line : nonFormatMetaData )
                mWriter.write(VCFHeader.METADATA_INDICATOR + line + "\n");

            // write out the column line
            StringBuilder b = new StringBuilder();
            b.append(VCFHeader.HEADER_INDICATOR);
            for (VCFHeader.HEADER_FIELDS field : header.getHeaderFields()) b.append(field + FIELD_SEPERATOR);
            if (header.hasGenotypingData()) {
                b.append("FORMAT" + FIELD_SEPERATOR);
                for (String field : header.getGenotypeSamples()) b.append(field + FIELD_SEPERATOR);
            }
            mWriter.write(b.toString() + "\n");
            mWriter.flush();  // necessary so that writing to an output stream will work
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
        addRecord(record, VCFGenotypeWriter.VALIDATION_STRINGENCY.STRICT);
    }

    /**
     * output a record to the VCF file
     *
     * @param record                the record to output
     * @param validationStringency  the validation stringency
     */
    public void addRecord(VCFRecord record, VCFGenotypeWriter.VALIDATION_STRINGENCY validationStringency) {
        if ( mHeader == null )
            throw new IllegalStateException("The VCF Header must be written before records can be added");

        String vcfString = record.toStringEncoding(mHeader, validationStringency);
        try {
            mWriter.write(vcfString + "\n");
            mWriter.flush();  // necessary so that writing to an output stream will work
        } catch (IOException e) {
            throw new RuntimeException("Unable to write the VCF object to a file");
        }

    }

    /**
     * attempt to close the VCF file
     */
    public void close() {
        try {
            mWriter.flush();
            mWriter.close();
        } catch (IOException e) {
            throw new RuntimeException("Unable to close VCFFile");
        }
    }

}
