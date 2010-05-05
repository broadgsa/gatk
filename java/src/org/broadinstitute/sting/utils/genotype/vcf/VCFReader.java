package org.broadinstitute.sting.utils.genotype.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import org.broad.tribble.FeatureReader;
import org.broad.tribble.index.linear.LinearIndex;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.refdata.tracks.builders.TribbleRMDTrackBuilder;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;

/** The VCFReader class, which given a valid vcf file, parses out the header and VCF records */
public class VCFReader implements Iterator<VCFRecord>, Iterable<VCFRecord> {

    // our VCF header
    private VCFHeader mHeader;

    // our iterator
    private Iterator<VCFRecord> iterator;

    // our feature reader; so we can close it
    private FeatureReader<VCFRecord> vcfReader = null;

    /**
     * Create a VCF reader, given a VCF file
     *
     * @param vcfFile the vcf file to write
     */
    public VCFReader(File vcfFile) {
        initialize(vcfFile, null);
    }

    /**
     * Create a VCF reader, given a VCF file
     *
     * @param vcfFile the vcf file to write
     */
    public VCFReader(File vcfFile, VCFCodec.LineTransform transform) {
        initialize(vcfFile, transform);
    }

    private void initialize(File vcfFile, VCFCodec.LineTransform transform) {
        VCFCodec codec = new VCFCodec();
        LinearIndex index = null;
        if (TribbleRMDTrackBuilder.requireIndex(vcfFile)) {
            try {
                index = TribbleRMDTrackBuilder.createIndex(vcfFile, new VCFCodec());
            } catch (IOException e) {
                throw new StingException("Unable to make required index for file " + vcfFile + " do you have write permissions to the directory?");
            }
        }
        if (transform != null) codec.setTransformer(transform);
        try {
            vcfReader = new FeatureReader(vcfFile,codec);
            iterator= vcfReader.iterator();
        } catch (FileNotFoundException e) {
            throw new StingException("Unable to read VCF File from " + vcfFile, e);
        } catch (IOException e) {
            throw new StingException("Unable to read VCF File from " + vcfFile, e);
        }
        mHeader = codec.getHeader();
    }


    /** @return true if we have another VCF record to return */
    public boolean hasNext() {
        return (iterator.hasNext());
    }

    /**
     * return the next available VCF record.  Make sure to check availability with a call to hasNext!
     *
     * @return a VCFRecord, representing the next record in the file
     */
    public VCFRecord next() {
        return iterator.next();
    }

    /** Remove is not supported */
    public void remove() {
        throw new UnsupportedOperationException("Unsupported operation");
    }



    /** @return get the header associated with this reader */
    public VCFHeader getHeader() {
        return this.mHeader;
    }

    public Iterator<VCFRecord> iterator() {
        return this;
    }

    public void close() {
        if (vcfReader != null) try {
            vcfReader.close();
        } catch (IOException e) {
            throw new StingException("Unable to close vcfReader",e);
        }
        iterator = null;
    }

}
