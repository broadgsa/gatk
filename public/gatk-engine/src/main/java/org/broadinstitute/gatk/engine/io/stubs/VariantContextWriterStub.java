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

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.index.IndexCreator;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.arguments.GATKArgumentCollection;
import org.broadinstitute.gatk.engine.io.OutputTracker;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import java.io.File;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.EnumSet;
import java.util.List;

/**
 * A stub for routing and management of genotype reading and writing.
 *
 * @author ebanks
 * @version 0.1
 */
public class VariantContextWriterStub implements Stub<VariantContextWriter>, VariantContextWriter {
    public final static boolean UPDATE_CONTIG_HEADERS = true;

    /**
     * The engine, central to the GATK's processing.
     */
    private final GenomeAnalysisEngine engine;

    /**
     * The file that this stub should write to.  Should be mutually
     * exclusive with genotypeStream.
     */
    private final File genotypeFile;

    /**
     * The output stream to which stub data should be written.  Will be
     * mutually exclusive with genotypeFile.
     */
    private final PrintStream genotypeStream;

    /**
     * A hack: push the argument sources into the VCF header so that the VCF header
     * can rebuild the command-line arguments.
     */
    private final Collection<Object> argumentSources;

    /**
     * Which IndexCreator to use
     */
    private final IndexCreator indexCreator;

    /**
     * The cached VCF header (initialized to null)
     */
    private VCFHeader vcfHeader = null;

    /**
     * Should we emit a compressed output stream?
     */
    private boolean isCompressed = false;

    /**
     * Should the header be written out?  A hidden argument.
     */
    private boolean skipWritingCommandLineHeader = false;

    /**
     * Should we not write genotypes even when provided?
     */
    private boolean doNotWriteGenotypes = false;

    /**
     * Should we force BCF writing regardless of the file extension?
     */
    private boolean forceBCF = false;

    /**
     * Should we write all of the fields in the FORMAT field, even if missing fields could be trimmed?
     */
    private boolean writeFullFormatField = false;

    /**
     * Connects this stub with an external stream capable of serving the
     * requests of the consumer of this stub.
     */
    protected OutputTracker outputTracker = null;

    /**
     * Create a new stub given the requested file.
     *
     * @param engine engine.
     * @param genotypeFile  file to (ultimately) create.
     * @param argumentSources sources.
     */
    public VariantContextWriterStub(GenomeAnalysisEngine engine, File genotypeFile, Collection<Object> argumentSources) {
        this.engine = engine;
        this.genotypeFile = genotypeFile;
        this.genotypeStream = null;

        this.indexCreator = GATKVCFUtils.makeIndexCreator(engine.getArguments().variant_index_type, engine.getArguments().variant_index_parameter,
                genotypeFile, null);
        this.argumentSources = argumentSources;
    }

    /**
     * Create a new stub given the requested file.
     *
     * @param engine engine.
     * @param genotypeStream  stream to (ultimately) write.
     * @param argumentSources sources.
     */
    public VariantContextWriterStub(GenomeAnalysisEngine engine, OutputStream genotypeStream, Collection<Object> argumentSources) {
        this.engine = engine;
        this.genotypeFile = null;
        this.genotypeStream = new PrintStream(genotypeStream);
        this.indexCreator = null;
        this.argumentSources = argumentSources;
    }

    /**
     * Retrieves the file to (ultimately) be created.
     * @return The file.  Can be null if genotypeStream is not.
     */
    public File getOutputFile() {
        return genotypeFile;
    }

    /**
     * Retrieves the output stream to which to (ultimately) write.
     * @return The file.  Can be null if genotypeFile is not.
     */
    public OutputStream getOutputStream() {
        return genotypeStream;
    }

    public boolean isCompressed() {
        return isCompressed;
    }

    public void setCompressed(final boolean compressed) {
        isCompressed = compressed;
    }

    public void setSkipWritingCommandLineHeader(final boolean skipWritingCommandLineHeader) {
        this.skipWritingCommandLineHeader = skipWritingCommandLineHeader;
    }

    public void setDoNotWriteGenotypes(final boolean doNotWriteGenotypes) {
        this.doNotWriteGenotypes = doNotWriteGenotypes;
    }

    public void setForceBCF(final boolean forceBCF) {
        this.forceBCF = forceBCF;
    }

    public void setWriteFullFormatField(final boolean writeFullFormatField) {
        this.writeFullFormatField = writeFullFormatField;
    }

    public IndexCreator getIndexCreator() {
        return indexCreator;
    }

    /**
     * Gets the master sequence dictionary from the engine associated with this stub
     * @link GenomeAnalysisEngine.getMasterSequenceDictionary
     * @return the master sequence dictionary from the engine associated with this stub
     */
    public SAMSequenceDictionary getMasterSequenceDictionary() {
        return engine.getMasterSequenceDictionary();
    }

    public EnumSet<Options> getWriterOptions() {
        return getWriterOptions(false);
    }

    public EnumSet<Options> getWriterOptions(boolean indexOnTheFly) {
        final List<Options> options = new ArrayList<>();

        if ( doNotWriteGenotypes ) options.add(Options.DO_NOT_WRITE_GENOTYPES);
        if ( engine.lenientVCFProcessing() ) options.add(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
        if ( indexOnTheFly) options.add(Options.INDEX_ON_THE_FLY);
        if ( writeFullFormatField ) options.add(Options.WRITE_FULL_FORMAT_FIELD);

        if ( forceBCF || (getOutputFile() != null && VariantContextWriterFactory.isBCFOutput(getOutputFile())) )
            options.add(Options.FORCE_BCF);

        return options.isEmpty() ? EnumSet.noneOf(Options.class) : EnumSet.copyOf(options);
    }

    /**
     * Retrieves the header to use when creating the new file.
     * @return header to use when creating the new file.
     */
    public VCFHeader getVCFHeader() {
        return vcfHeader;
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
        setDoNotWriteGenotypes(argumentCollection.sitesOnlyVCF);
        setSkipWritingCommandLineHeader(argumentCollection.disableCommandLineInVCF);
        setForceBCF(argumentCollection.forceBCFOutput);
        setWriteFullFormatField(argumentCollection.neverTrimVCFFormatField);
    }

    @Override
    public void writeHeader(VCFHeader header) {
        vcfHeader = header;

        if ( header.isWriteEngineHeaders() ) {
            // skip writing the command line header if requested
            if ( ! skipWritingCommandLineHeader && header.isWriteCommandLine() ) {
                // Always add the header line, as the current format allows multiple entries
                final VCFHeaderLine commandLineArgHeaderLine = GATKVCFUtils.getCommandLineArgumentHeaderLine(vcfHeader, engine, argumentSources);
                vcfHeader.addMetaDataLine(commandLineArgHeaderLine);
            }

            if ( UPDATE_CONTIG_HEADERS )
                vcfHeader = GATKVCFUtils.withUpdatedContigs(vcfHeader, engine);
        }

        outputTracker.getStorage(this).writeHeader(vcfHeader);
    }

    /**
     * @{inheritDoc}
     */
    public void add(VariantContext vc) {
        outputTracker.getStorage(this).add(vc);
    }

    /**
     * @{inheritDoc}
     */
    @Override
    public void close() {
        outputTracker.getStorage(this).close();
    }

    /**
     * Gets a string representation of this object.
     * @return a string representation of this object.
     */
    @Override
    public String toString() {
        return (getOutputFile() == null) ? "(Stream)" : getOutputFile().getAbsolutePath();
    }

    /**
     * Should we also write a BCF file alongside our VCF file for testing
     *
     * TODO -- remove me when argument generateShadowBCF is removed
     *
     * @return
     */
    public boolean alsoWriteBCFForTest() {
        return engine.getArguments().numberOfDataThreads == 1 && // only works single threaded
                ! isCompressed() && // for non-compressed outputs
                getOutputFile() != null && // that are going to disk
                engine.getArguments().generateShadowBCF; // and we actually want to do it
    }

    /**
     * Check the return from PrintStream.checkError() if underlying stream for a java.io.PrintStream
     * @return true if PrintStream.checkError() returned true, false otherwise
     */
    public boolean checkError(){
        return genotypeStream.checkError();
    }
}
