/*
 * Copyright (c) 2012, The Broad Institute
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

import net.sf.samtools.SAMSequenceDictionary;
import org.broadinstitute.sting.gatk.CommandLineExecutable;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.io.OutputTracker;
import org.broadinstitute.sting.utils.classloader.JVMUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.variantcontext.writer.Options;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

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
     * The cached VCF header (initialized to null)
     */
    private VCFHeader vcfHeader = null;

    /**
     * Should we emit a compressed output stream?
     */
    private final boolean isCompressed;

    /**
     * A hack: push the argument sources into the VCF header so that the VCF header
     * can rebuild the command-line arguments.
     */
    private final Collection<Object> argumentSources;

    /**
     * Should the header be written out?  A hidden argument.
     */
    private final boolean skipWritingCommandLineHeader;

    /**
     * Should we not write genotypes even when provided?
     */
    private final boolean doNotWriteGenotypes;

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
     * @param isCompressed  should we compress the output stream?
     * @param argumentSources sources.
     * @param skipWritingCommandLineHeader skip writing header.
     * @param doNotWriteGenotypes do not write genotypes.
     */
    public VariantContextWriterStub(GenomeAnalysisEngine engine, File genotypeFile, boolean isCompressed, Collection<Object> argumentSources, boolean skipWritingCommandLineHeader, boolean doNotWriteGenotypes) {
        this.engine = engine;
        this.genotypeFile = genotypeFile;
        this.genotypeStream = null;
        this.isCompressed = isCompressed;
        this.argumentSources = argumentSources;
        this.skipWritingCommandLineHeader = skipWritingCommandLineHeader;
        this.doNotWriteGenotypes = doNotWriteGenotypes;
    }

    /**
     * Create a new stub given the requested file.
     *
     * @param engine engine.
     * @param genotypeStream  stream to (ultimately) write.
     * @param isCompressed  should we compress the output stream?
     * @param argumentSources sources.
     * @param skipWritingCommandLineHeader skip writing header.
     * @param doNotWriteGenotypes do not write genotypes.
     */
    public VariantContextWriterStub(GenomeAnalysisEngine engine, OutputStream genotypeStream, boolean isCompressed, Collection<Object> argumentSources, boolean skipWritingCommandLineHeader, boolean doNotWriteGenotypes) {
        this.engine = engine;
        this.genotypeFile = null;
        this.genotypeStream = new PrintStream(genotypeStream);
        this.isCompressed = isCompressed;
        this.argumentSources = argumentSources;
        this.skipWritingCommandLineHeader = skipWritingCommandLineHeader;
        this.doNotWriteGenotypes = doNotWriteGenotypes;
    }

    /**
     * Retrieves the file to (ultimately) be created.
     * @return The file.  Can be null if genotypeStream is not.
     */
    public File getFile() {
        return genotypeFile;
    }

    /**
     * Retrieves the output stearm to which to (ultimately) write.
     * @return The file.  Can be null if genotypeFile is not.
     */
    public PrintStream getOutputStream() {
        return genotypeStream;
    }

    /**
     * Retrieves the output stearm to which to (ultimately) write.
     * @return The file.  Can be null if genotypeFile is not.
     */
    public boolean isCompressed() {
        return isCompressed;
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
        List<Options> options = new ArrayList<Options>();

        if ( doNotWriteGenotypes ) options.add(Options.DO_NOT_WRITE_GENOTYPES);
        if ( engine.getArguments().allowMissingVCFHeaders ) options.add(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
        if ( indexOnTheFly && ! isCompressed() ) options.add(Options.INDEX_ON_THE_FLY);

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

    public void writeHeader(VCFHeader header) {
        vcfHeader = header;

        if ( header.isWriteEngineHeaders() ) {
            // skip writing the command line header if requested
            if ( ! skipWritingCommandLineHeader && header.isWriteCommandLine() ) {
                // Check for the command-line argument header line.  If not present, add it in.
                final VCFHeaderLine commandLineArgHeaderLine = getCommandLineArgumentHeaderLine();
                final boolean foundCommandLineHeaderLine = vcfHeader.getMetaDataLine(commandLineArgHeaderLine.getKey()) != null;
                if ( ! foundCommandLineHeaderLine )
                    vcfHeader.addMetaDataLine(commandLineArgHeaderLine);
            }

            if ( UPDATE_CONTIG_HEADERS )
                vcfHeader = VCFUtils.withUpdatedContigs(vcfHeader, engine);
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
    public void close() {
        outputTracker.getStorage(this).close();
    }

    /**
     * Gets a string representation of this object.
     * @return a string representation of this object.
     */
    @Override
    public String toString() {
        return getClass().getName();
    }

    /**
     * Should we also write a BCF file alongside our VCF file for testing
     *
     * TODO -- remove me when argument generateShadowBCF is removed
     *
     * @return
     */
    public boolean alsoWriteBCFForTest() {
        return engine.getArguments().numberOfThreads == 1 && // only works single threaded
                ! isCompressed() && // for non-compressed outputs
                getFile() != null && // that are going to disk
                engine.getArguments().generateShadowBCF; // and we actually want to do it
    }

    /**
     * Gets the appropriately formatted header for a VCF file
     * @return VCF file header.
     */
    private VCFHeaderLine getCommandLineArgumentHeaderLine() {
        CommandLineExecutable executable = JVMUtils.getObjectOfType(argumentSources,CommandLineExecutable.class);
        return new VCFHeaderLine(executable.getAnalysisName(), "\"" + engine.createApproximateCommandLineArgumentString(argumentSources.toArray()) + "\"");
    }
}
