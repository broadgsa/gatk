/*
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.io.stubs;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileHeader;

import java.io.File;
import java.io.OutputStream;

import org.broadinstitute.sting.gatk.io.OutputTracker;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.StingException;

/**
 * A stub for routing and management of SAM file reading and writing.
 *
 * @author mhanna
 * @version 0.1
 */
public class SAMFileWriterStub implements Stub<SAMFileWriter>, StingSAMFileWriter {
    /**
     * Engine to use for collecting attributes for the output SAM file.
     */
    private final GenomeAnalysisEngine engine;

    /**
     * A header supplied by the user that overrides the merged header from the input BAM.
     */
    private SAMFileHeader headerOverride = null;

    /**
     * The sam file that this stub should write to.  Should be passed along to
     * whatever happens to create the StreamConnector.
     */
    private final File samFile;

    /**
     * The target output stream, to be used in place of the SAM file.
     */
    private final OutputStream samOutputStream;

    /**
     * The validation stringency to apply when reading this file.
     */
    private Integer compressionLevel = null;

    /**
     * Should this BAM be presorted?
     */
    private boolean presorted = true;

    /**
     * Connects this stub with an external stream capable of serving the
     * requests of the consumer of this stub.
     */
    private OutputTracker outputTracker = null;

    /**
     * Has the write started?  If so, throw an exception if someone tries to
     * change write parameters to the file (compression level, presorted flag,
     * header, etc).
     */
    private boolean writeStarted = false;

    /**
     * Create a new stub given the requested SAM file and compression level.
     * @param engine source of header data, maybe other data about input files.
     * @param samFile SAM file to (ultimately) create.
     */
    public SAMFileWriterStub( GenomeAnalysisEngine engine, File samFile ) {
        this.engine = engine;
        this.samFile = samFile;
        this.samOutputStream = null;
    }

    /**
     * Create a new stub given the requested SAM file and compression level.
     * @param engine source of header data, maybe other data about input files.
     * @param stream Output stream to which data should be written.
     */
    public SAMFileWriterStub( GenomeAnalysisEngine engine, OutputStream stream ) {
        this.engine = engine;
        this.samFile = null;
        this.samOutputStream = stream;
    }

    /**
     * Retrieves the SAM file to (ultimately) be created.
     * @return The SAM file.  Must not be null.
     */
    public File getSAMFile() {
        return samFile;
    }

    public OutputStream getSAMOutputStream() {
        return samOutputStream;
    }

    /**
     * Retrieves the header to use when creating the new SAM file.
     * @return header to use when creating the new SAM file.
     */
    public SAMFileHeader getFileHeader() {
        return headerOverride != null ? headerOverride : engine.getSAMFileHeader();
    }

    /**
     * Retrieves the desired compression level for 
     * @return The current compression level.  Could be null if the user doesn't care.
     */
    public Integer getCompressionLevel() {
        return compressionLevel;
    }

    /**
     * Sets the desired compression level.
     * @param compressionLevel The suggested compression level.
     */
    public void setCompressionLevel( Integer compressionLevel ) {
        if(writeStarted)
            throw new StingException("User attempted to change the compression level of a file with alignments already in it.");        
        this.compressionLevel = compressionLevel;
    }

    /**
     * Whether the BAM file to create is actually presorted.
     * @return True if the BAM file is presorted.  False otherwise.
     */
    public boolean isPresorted() {
        return this.presorted;
    }

    /**
     * Set Whether the BAM file to create is actually presorted.
     * @param presorted True if the BAM file is presorted.  False otherwise.
     */
    public void setPresorted(boolean presorted) {
        if(writeStarted)
            throw new StingException("User attempted to change the presorted state of a file with alignments already in it.");
        this.presorted = presorted;
    }

    /**
     * Registers the given streamConnector with this stub.
     * @param outputTracker The connector used to provide an appropriate stream.
     */
    public void register( OutputTracker outputTracker ) {
        this.outputTracker = outputTracker;
    }

    /**
     * Use the given header as the target for this writer.
     * @param header The header to write.
     */
    public void writeHeader(SAMFileHeader header) {
        if(writeStarted)
            throw new StingException("User attempted to change the header of a file with alignments already in it.");
        this.headerOverride = header;
    }

    /**
     * @{inheritDoc}
     */
    public void addAlignment( SAMRecord alignment ) {
        writeStarted = true;
        outputTracker.getStorage(this).addAlignment(alignment);
    }

    /**
     * @{inheritDoc}
     */
    public void close() {
        outputTracker.getStorage(this).close();    
    }

}
