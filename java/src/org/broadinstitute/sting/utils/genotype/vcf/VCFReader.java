package org.broadinstitute.sting.utils.genotype.vcf;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;

import org.broad.tribble.FeatureReader;
import org.broad.tribble.index.Index;
import org.broad.tribble.readers.BasicFeatureReader;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.refdata.tracks.builders.TribbleRMDTrackBuilder;
import org.broadinstitute.sting.utils.StingException;

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
        initialize(vcfFile, null, true);
    }

    /**
     * Create a VCF reader, given a VCF file
     *
     * @param vcfFile the vcf file to write
     */
    public VCFReader(File vcfFile, boolean createIndexOnDisk) {
        initialize(vcfFile, null, createIndexOnDisk);
    }

    /**
     * Create a VCF reader, given a VCF file
     *
     * @param vcfFile the vcf file to write
     */
    public VCFReader(File vcfFile, VCFCodec.LineTransform transform) {
        initialize(vcfFile, transform, true);
    }

    /**
     * initialize the VCF reader
     * @param vcfFile the VCF file to open
     * @param transform the line transformer to use, if any
     * @param createIndexOnDisk do we need to create an index on disk?
     */
    private void initialize(File vcfFile, VCFCodec.LineTransform transform, boolean createIndexOnDisk) {
        VCFCodec codec = new VCFCodec();
        Index index = createIndex(vcfFile, createIndexOnDisk);
        if (transform != null) codec.setTransformer(transform);
        try {
            vcfReader = new BasicFeatureReader(vcfFile.getAbsolutePath(),index,codec);
            iterator= vcfReader.iterator();
        } catch (FileNotFoundException e) {
            throw new StingException("Unable to read VCF File from " + vcfFile, e);
        } catch (IOException e) {
            throw new StingException("Unable to read VCF File from " + vcfFile, e);
        }
        mHeader = codec.getHeader();
    }

    /**
     * create an index given:
     * @param vcfFile the vcf file
     * @param createIndexOnDisk do we create the index on disk (or only in memory?)
     * @return an instance of an index
     */
    private Index createIndex(File vcfFile, boolean createIndexOnDisk) {
        Index index = null;
        try {
            index = TribbleRMDTrackBuilder.loadIndex(vcfFile, new VCFCodec(), createIndexOnDisk);
        } catch (IOException e) {
            throw new StingException("Unable to make required index for file " + vcfFile + " do you have write permissions to the directory?");
        }

        return index;
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

    /**
     * close the files
     */
    public void close() {
        if (vcfReader != null) try {
            vcfReader.close();
        } catch (IOException e) {
            throw new StingException("Unable to close vcfReader",e);
        }
        iterator = null;
    }

}
