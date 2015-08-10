/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.io.stubs;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.ProgressLoggerInterface;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.arguments.GATKArgumentCollection;
import org.broadinstitute.gatk.engine.io.OutputTracker;
import org.broadinstitute.gatk.utils.io.ReferenceBacked;
import org.broadinstitute.gatk.utils.sam.GATKSAMFileWriter;
import org.broadinstitute.gatk.engine.iterators.ReadTransformer;
import org.broadinstitute.gatk.utils.baq.BAQ;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.io.File;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;

/**
 * A stub for routing and management of SAM file reading and writing.
 *
 * @author mhanna
 * @version 0.1
 */
public class SAMFileWriterStub implements Stub<SAMFileWriter>, GATKSAMFileWriter, ReferenceBacked {
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
     * The reference file for stub.
     */
    private File referenceFile;

    /**
     * The target output stream, to be used in place of the SAM file.
     */
    private final OutputStream samOutputStream;

    /**
     * The validation stringency to apply when reading this file.
     */
    private Integer compressionLevel = null;

    /**
     * Should the GATK index the output BAM on-the-fly?
     */
    private boolean indexOnTheFly = false;

    /**
     * Should the GATK generate an md5 for the output BAM?
     */
    private boolean generateMD5 = false;

    /**
     * Should this BAM be presorted?
     */
    private boolean presorted = true;

    /**
     * How many records should the BAM writer store in RAM while
     * sorting the BAM on-the-fly?
     */
    private Integer maxRecordsInRam = null;

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
     * HMM for BAQ, if needed
     */
    BAQ baqHMM = new BAQ();

    /**
     * Should we simplify the BAM file while writing it out?
     */
    private boolean simplifyBAM = false;

    private List<ReadTransformer> onOutputReadTransformers = null;

    /**
     * Create a new stub given the requested SAM file and compression level.
     * @param engine source of header data, maybe other data about input files.
     * @param samFile SAM file to (ultimately) create.
     */
    public SAMFileWriterStub( GenomeAnalysisEngine engine, File samFile ) {
        this(engine, samFile, null);
    }

    /**
     * Create a new stub given the requested SAM file and compression level.
     * @param engine source of header data, maybe other data about input files.
     * @param stream Output stream to which data should be written.
     */
    public SAMFileWriterStub( GenomeAnalysisEngine engine, OutputStream stream ) {
        this(engine, null, stream);
    }

    private SAMFileWriterStub(final GenomeAnalysisEngine engine, final File samFile, final OutputStream stream) {
        this.engine = engine;
        this.samFile = samFile;
        this.samOutputStream = stream;
    }

    /**
     * Creates a SAMFileWriter using all of the features currently set in the engine (command line arguments, ReadTransformers, etc)
     * @param file the filename to write to
     * @param engine the engine
     * @return a SAMFileWriter with the correct options set
     */
    public static SAMFileWriter createSAMFileWriter(final String file, final GenomeAnalysisEngine engine) {
        final SAMFileWriterStub output = new SAMFileWriterStub(engine, new File(file));
        output.processArguments(engine.getArguments());
        return output;
    }

    /**
     *  As {@link #createSAMFileWriter(String, org.broadinstitute.gatk.engine.GenomeAnalysisEngine)}, but also sets the header
     */
    public static SAMFileWriter createSAMFileWriter(final String file, final GenomeAnalysisEngine engine, final SAMFileHeader header) {
        final SAMFileWriterStub output = (SAMFileWriterStub) createSAMFileWriter(file, engine);
        output.writeHeader(header);
        return output;
    }

    /**
     * Retrieves the SAM file to (ultimately) be created.
     * @return The SAM file.  Must not be null.
     */
    public File getOutputFile() {
        return samFile;
    }

    public boolean simplifyBAM() {
        return simplifyBAM;
    }

    public void setSimplifyBAM(boolean v) {
        simplifyBAM = v;
    }

    public OutputStream getOutputStream() {
        return samOutputStream;
    }

    @Override
    public File getReferenceFile() {
        return referenceFile;
    }

    @Override
    public void setReferenceFile(final File referenceFile) {
        this.referenceFile = referenceFile;
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
            throw new ReviewedGATKException("Attempted to change the compression level of a file with alignments already in it.");
        this.compressionLevel = compressionLevel;
    }

    /**
     * Gets whether to index this output stream on-the-fly.
     * @return True means create an index.  False means skip index creation.
     */
    public Boolean getIndexOnTheFly() {
        return indexOnTheFly;
    }

    /**
     * Controls whether to index this output stream on-the-fly.
     * @param indexOnTheFly True means create an index.  False means skip index creation.
     */
    public void setIndexOnTheFly( boolean indexOnTheFly ) {
        if(writeStarted)
            throw new UserException("Attempted to index a BAM on the fly of a file with alignments already in it.");
        this.indexOnTheFly = indexOnTheFly;
    }

    /**
     * Gets whether to generate an md5 on-the-fly for this BAM.
     * @return True generates the md5.  False means skip writing the file.
     */
    public Boolean getGenerateMD5() {
        return generateMD5;
    }

    /**
     * Gets whether to generate an md5 on-the-fly for this BAM.
     * @param generateMD5   True generates the md5.  False means skip writing the file.
     */
    public void setGenerateMD5(boolean generateMD5) {
        if(writeStarted)
            throw new UserException("Attempted to turn on md5 generation for BAM file with alignments already in it.");        
        this.generateMD5 = generateMD5;
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
            throw new ReviewedGATKException("Attempted to change the presorted state of a file with alignments already in it.");
        this.presorted = presorted;
    }

    /**
     * Get the maximum number of reads to hold in RAM when sorting a BAM on-the-fly.
     * @return Max records in RAM, or null if unset.
     */
    public Integer getMaxRecordsInRam() {
        return this.maxRecordsInRam;
    }

    /**
     * Sets the maximum number of reads to hold in RAM when sorting a BAM on-the-fly.
     * @param maxRecordsInRam Max number of records in RAM.
     */
    public void setMaxRecordsInRam(int maxRecordsInRam) {
        if(writeStarted)
            throw new ReviewedGATKException("Attempted to change the max records in RAM of a file with alignments already in it.");
        this.maxRecordsInRam = maxRecordsInRam;
    }

    /**
     * Registers the given streamConnector with this stub.
     * @param outputTracker The connector used to provide an appropriate stream.
     */
    public void register( OutputTracker outputTracker ) {
        this.outputTracker = outputTracker;
    }

    @Override
    public void processArguments( final GATKArgumentCollection argumentCollection ) {
        if (argumentCollection.bamCompression != null)
            setCompressionLevel(argumentCollection.bamCompression);
        setGenerateMD5(argumentCollection.enableBAMmd5);
        setIndexOnTheFly(!argumentCollection.disableBAMIndexing);
        setSimplifyBAM(argumentCollection.simplifyBAM);

    }

    /**
     * Use the given header as the target for this writer.
     * @param header The header to write.
     */
    public void writeHeader(SAMFileHeader header) {
        if(writeStarted)
            throw new ReviewedGATKException("Attempted to change the header of a file with alignments already in it.");
        this.headerOverride = header;
    }

    private void initializeReadTransformers() {
        this.onOutputReadTransformers = new ArrayList<>(engine.getReadTransformers().size());
        for ( final ReadTransformer transformer : engine.getReadTransformers() ) {
            if ( transformer.getApplicationTime() == ReadTransformer.ApplicationTime.ON_OUTPUT )
                onOutputReadTransformers.add(transformer);
        }
    }

    /**
     * @{inheritDoc}
     */
    public void addAlignment( final SAMRecord readIn ) {
        if ( onOutputReadTransformers == null )
            initializeReadTransformers();

        GATKSAMRecord workingRead = (GATKSAMRecord)readIn;

        // run on output read transformers
        for ( final ReadTransformer transform : onOutputReadTransformers )
            workingRead = transform.apply(workingRead);

        writeStarted = true;
        outputTracker.getStorage(this).addAlignment(workingRead);
    }

    /**
     * @{inheritDoc}
     */
    public void close() {
        outputTracker.getStorage(this).close();    
    }

    /**
     * @throws java.lang.UnsupportedOperationException No progress logging in this implementation.
     */
    @Override
    public void setProgressLogger(final ProgressLoggerInterface logger) {
        throw new UnsupportedOperationException("Progress logging not supported");
    }
}
