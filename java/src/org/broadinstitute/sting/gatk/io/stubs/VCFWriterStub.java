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

import java.io.File;
import java.io.PrintStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.TreeSet;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFWriter;
import org.broadinstitute.sting.commandline.CommandLineUtils;
import org.broadinstitute.sting.gatk.CommandLineExecutable;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.io.OutputTracker;
import org.broadinstitute.sting.utils.classloader.JVMUtils;

/**
 * A stub for routing and management of genotype reading and writing.
 *
 * @author ebanks
 * @version 0.1
 */
public class VCFWriterStub implements Stub<VCFWriter>, VCFWriter {
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
    private final boolean skipWritingHeader;

    /**
     * Connects this stub with an external stream capable of serving the
     * requests of the consumer of this stub.
     */
    protected OutputTracker outputTracker = null;

    /**
     * Create a new stub given the requested file.
     * @param genotypeFile  file to (ultimately) create.
     * @param isCompressed  should we compress the output stream?
     */
    public VCFWriterStub(GenomeAnalysisEngine engine, File genotypeFile, boolean isCompressed, Collection<Object> argumentSources, boolean skipWritingHeader) {
        this.engine = engine;
        this.genotypeFile = genotypeFile;
        this.genotypeStream = null;
        this.isCompressed = isCompressed;
        this.argumentSources = argumentSources;
        this.skipWritingHeader = skipWritingHeader;
    }

    /**
     * Create a new stub given the requested file.
     * @param genotypeStream  stream to (ultimately) write.
     * @param isCompressed  should we compress the output stream?
     */
    public VCFWriterStub(GenomeAnalysisEngine engine, OutputStream genotypeStream, boolean isCompressed, Collection<Object> argumentSources, boolean skipWritingHeader) {
        this.engine = engine;
        this.genotypeFile = null;
        this.genotypeStream = new PrintStream(genotypeStream);
        this.isCompressed = isCompressed;
        this.argumentSources = argumentSources;
        this.skipWritingHeader = skipWritingHeader;
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

        // Check for the command-line argument header line.  If not present, add it in.
        VCFHeaderLine commandLineArgHeaderLine = getCommandLineArgumentHeaderLine();
        boolean foundCommandLineHeaderLine = false;
        for(VCFHeaderLine line: vcfHeader.getMetaData()) {
            if(line.getKey().equals(commandLineArgHeaderLine.getKey()))
                foundCommandLineHeaderLine = true;
        }
        if(!foundCommandLineHeaderLine && !skipWritingHeader)
            vcfHeader.addMetaDataLine(commandLineArgHeaderLine);

        outputTracker.getStorage(this).writeHeader(vcfHeader);
    }

    /**
     * @{inheritDoc}
     */
    public void add(VariantContext vc, byte ref) {
        outputTracker.getStorage(this).add(vc,ref);
    }

    /**
     * @{inheritDoc}
     */
    public void close() {
        outputTracker.getStorage(this).close();
    }

    /**
     * Gets a string representation of this object.
     * @return
     */
    @Override
    public String toString() {
        return getClass().getName();
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