/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package net.sf.samtools;

import net.sf.samtools.util.CloseableIterator;
import net.sf.picard.PicardException;

import java.io.*;
import java.util.List;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.lang.reflect.InvocationTargetException;

import org.broadinstitute.sting.utils.JVMUtils;
import org.broadinstitute.sting.utils.StingException;

/**
 * Class for reading and querying SAM/BAM files.  Delegates to appropriate concrete implementation.
 */
public class SAMFileReader2 extends SAMFileReader {
    /**
     * Prepare to read a SAM or BAM file.  If the given file is a BAM, and has a companion BAI index file
     * that is named according to the convention, it will be found and opened, and indexed query will be allowed.
     */
    public SAMFileReader2(final File file) {
        this(file, null, false);
    }

    /**
     * Read a SAM or BAM file, possibly with an index file if present.
     * If the given file is a BAM, and an index is present, indexed query will be allowed.
     *
     * @param file SAM or BAM.
     * @param eagerDecode if true, decode SAM record entirely when reading it.
     */
    public SAMFileReader2(final File file, final boolean eagerDecode) {
        this(file,null,eagerDecode);
    }

    /**
     * Read a SAM or BAM file, possibly with an index file. If the given file is a BAM, and an index is present,
     * indexed query will be allowed.
     *
     * @param file SAM or BAM.
     * @param indexFile Location of index file, or null in order to use the default index file (if present).
     * @param eagerDecode eagerDecode if true, decode SAM record entirely when reading it.
     */
    public SAMFileReader2(final File file, final File indexFile, final boolean eagerDecode){
        super(file,indexFile,eagerDecode);
        close();

        try {
            BAMFileReader2 reader = new BAMFileReader2(file,eagerDecode,getDefaultValidationStringency());
            BAMFileIndex2 index = new BAMFileIndex2(indexFile != null ? indexFile : findIndexFileFromParent(file));
            reader.setFileIndex(index);

            JVMUtils.setFieldValue(getField("mReader"),this,reader);
            JVMUtils.setFieldValue(getField("mFileIndex"),this,index);
        }
        catch(IOException ex) {
            throw new StingException("Unable to load BAM file: " + file,ex);
        }
    }

    /**
     * Get the number of levels employed by this index.
     * @return Number of levels in this index.
     */
    public int getNumIndexLevels() {
        final BAMFileIndex2 fileIndex = (BAMFileIndex2)JVMUtils.getFieldValue(getField("mFileIndex"),this);
        if(fileIndex == null)
            throw new SAMException("Unable to determine number of index levels; BAM file index is not present.");
        return fileIndex.getNumIndexLevels();
    }

    /**
     * Gets the level associated with the given bin number.
     * @param bin The bin for which to determine the level.
     * @return the level associated with the given bin number.
     */
    public int getLevelForBin(final Bin bin) {
        final BAMFileIndex2 fileIndex = (BAMFileIndex2)JVMUtils.getFieldValue(getField("mFileIndex"),this);
        if(fileIndex == null)
            throw new SAMException("Unable to determine number of index levels; BAM file index is not present.");
        return fileIndex.getLevelForBinNumber(bin.binNumber);
    }

    /**
     * Gets the first locus that this bin can index into.
     * @param bin The bin to test.
     * @return The last position that the given bin can represent.
     */
    public int getFirstLocusInBin(final Bin bin) {
        final BAMFileIndex2 fileIndex = (BAMFileIndex2)JVMUtils.getFieldValue(getField("mFileIndex"),this);
        if(fileIndex == null)
            throw new SAMException("Unable to determine number of index levels; BAM file index is not present.");
        return fileIndex.getFirstLocusInBin(bin);
    }

    /**
     * Gets the last locus that this bin can index into.
     * @param bin The bin to test.
     * @return The last position that the given bin can represent.
     */
    public int getLastLocusInBin(final Bin bin) {
        final BAMFileIndex2 fileIndex = (BAMFileIndex2)JVMUtils.getFieldValue(getField("mFileIndex"),this);
        if(fileIndex == null)
            throw new SAMException("Unable to determine number of index levels; BAM file index is not present.");
        return fileIndex.getLastLocusInBin(bin);
    }

    /**
     * Iterate through the given chunks in the file.
     * @param chunks List of chunks for which to retrieve data.
     * @return An iterator over the given chunks.
     */
    public CloseableIterator<SAMRecord> iterator(List<Chunk> chunks) {
        // TODO: Add sanity checks so that we're not doing this against an unsupported BAM file.
        BAMFileReader2 reader = (BAMFileReader2)JVMUtils.getFieldValue(getField("mReader"),this);
        return reader.getIterator(chunks);
    }

    public List<Bin> getOverlappingBins(final String sequence, final int start, final int end) {
        // TODO: Add sanity checks so that we're not doing this against an unsupported BAM file.
        BAMFileReader2 reader = (BAMFileReader2)JVMUtils.getFieldValue(getField("mReader"),this);
        return reader.getOverlappingBins(sequence,start,end);
    }

    public List<Chunk> getFilePointersBounding(final String sequence, final int start, final int end) {
        // TODO: Add sanity checks so that we're not doing this against an unsupported BAM file.
        BAMFileReader2 reader = (BAMFileReader2)JVMUtils.getFieldValue(getField("mReader"),this);
        return reader.getFilePointersBounding(sequence,start,end);
    }

    private Field getField(String fieldName) {
        try {
            return getClass().getSuperclass().getDeclaredField(fieldName);
        }
        catch(NoSuchFieldException ex) {
            throw new StingException("Unable to load field: " + fieldName);
        }
    }

    private File findIndexFileFromParent(File bamFile) {
        try {
            Method method = getClass().getSuperclass().getDeclaredMethod("findIndexFile",File.class);
            method.setAccessible(true);
            return (File)method.invoke(this,bamFile);
        }
        catch(IllegalAccessException ex) {
            throw new StingException("Unable to run method findIndexFile",ex);
        }
        catch(InvocationTargetException ex) {
            throw new StingException("Unable to run method findIndexFile",ex);
        }
        catch(NoSuchMethodException ex) {
            throw new StingException("Unable to run method findIndexFile",ex);
        }

    }
}
