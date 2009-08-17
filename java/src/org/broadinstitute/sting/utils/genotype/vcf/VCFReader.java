package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.utils.StingException;

import java.io.*;
import java.util.List;
import java.util.Iterator;
import java.util.ArrayList;
import java.nio.ByteBuffer;
import java.nio.charset.Charset;

/**
 * The VCFReader class, which given a valid vcf file, parses out the header and VCF records
 */
public class VCFReader implements Iterator<VCFRecord> {

    // our VCF header
    private VCFHeader mHeader;

    // our buffered input stream
    private BufferedReader mReader;

    // our next record
    private VCFRecord mNextRecord = null;

    /**
     * Create a VCF reader, given a VCF file
     *
     * @param vcfFile
     */
    public VCFReader(File vcfFile) {
        Charset utf8 = Charset.forName("UTF-8");
        try {
            mReader = new BufferedReader(
                    new InputStreamReader(
                            new FileInputStream(vcfFile),
                            utf8));
        } catch (FileNotFoundException e) {
            throw new StingException("Unable to find VCF file: " + vcfFile, e);
        }

        String line = null;
        try {
            ArrayList<String> lines = new ArrayList<String>();
            line = mReader.readLine();
            while (line.startsWith("#")) {
                lines.add(line);
                line = mReader.readLine();
            }
            mHeader = new VCFHeader(lines);
            mNextRecord = new VCFRecord(mHeader, line);
        } catch (IOException e) {
            throw new StingException("Failed to parse VCF File on line: " + line, e);
        }

    }

    /**
     *
     * @return true if we have another VCF record to return
     */
    public boolean hasNext() {
        return (mNextRecord != null);
    }

    /**
     * return the next available VCF record.  Make sure to check availability with a call to hasNext!
     * @return a VCFRecord, representing the next record in the file
     */
    public VCFRecord next() {
        VCFRecord rec = mNextRecord;
        try {
            String line = mReader.readLine();
            if (line == null) mNextRecord = null;
            else mNextRecord = new VCFRecord(mHeader, line);
        } catch (IOException e) {
            mNextRecord = null;
        }
        return rec;
    }

    /**
     * Remove is not supported
     */
    public void remove() {
        throw new UnsupportedOperationException("Unsupported operation");
    }
}
