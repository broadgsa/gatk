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

package org.broadinstitute.gatk.engine.io.storage;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import org.apache.log4j.Logger;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import org.broadinstitute.gatk.engine.io.stubs.VariantContextWriterStub;
import org.broadinstitute.gatk.utils.refdata.tracks.FeatureManager;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import htsjdk.variant.bcf2.BCF2Utils;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFHeader;

import java.io.*;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;

/**
 * Provides temporary and permanent storage for genotypes in VCF format.
 *
 * @author mhanna
 * @version 0.1
 */
public class VariantContextWriterStorage implements Storage<VariantContextWriterStorage>, VariantContextWriter {
    /**
     * our log, which we want to capture anything from this class
     */
    private static Logger logger = Logger.getLogger(VariantContextWriterStorage.class);

    private final static int BUFFER_SIZE = 1048576;

    protected final File file;
    protected OutputStream stream;
    protected final VariantContextWriter writer;
    boolean closed = false;

    /**
     * Constructs an object which will write directly into the output file provided by the stub.
     * Intentionally delaying the writing of the header -- this should be filled in by the walker.
     *
     * Respecs the isCompressed() request in stub, so if isCompressed() is true then this
     * will create a storage output that dumps output to a BlockCompressedOutputStream.
     *
     * @param stub Stub to use when constructing the output file.
     */
    public VariantContextWriterStorage(VariantContextWriterStub stub)  {
        if ( stub.getOutputFile() != null ) {
            this.file = stub.getOutputFile();
            writer = vcfWriterToFile(stub,stub.getOutputFile(),true,true);
        }
        else if ( stub.getOutputStream() != null ) {
            this.file = null;
            this.stream = stub.getOutputStream();
            writer = VariantContextWriterFactory.create(stream,
                    stub.getMasterSequenceDictionary(), stub.getWriterOptions(false));
        }
        else
            throw new ReviewedGATKException("Unable to create target to which to write; storage was provided with neither a file nor a stream.");
    }

    /**
     * Constructs an object which will redirect into a different file.
     *
     * Note that this function does not respect the isCompressed() request from the stub, in order
     * to ensure that tmp. files can be read back in by the Tribble system, and merged with the mergeInto function.
     *
     * @param stub Stub to use when synthesizing file / header info.
     * @param tempFile File into which to direct the output data.
     */
    public VariantContextWriterStorage(VariantContextWriterStub stub, File tempFile) {
        //logger.debug("Creating temporary output file " + tempFile.getAbsolutePath() + " for VariantContext output.");
        this.file = tempFile;
        this.writer = vcfWriterToFile(stub, file, false, false);
        writer.writeHeader(stub.getVCFHeader());
    }

    /**
     * common initialization routine for multiple constructors
     * @param stub Stub to use when constructing the output file.
     * @param file Target file into which to write VCF records.
     * @param indexOnTheFly true to index the file on the fly.  NOTE: will be forced to false for compressed files.
     * @param allowCompressed if false, we won't compress the output, even if the stub requests it.  Critical
     *                        for creating temp. output files that will be subsequently merged, as these do not
     *                        support compressed output
     * @return A VCF writer for use with this class
     */
    private VariantContextWriter vcfWriterToFile(final VariantContextWriterStub stub,
                                                 final File file,
                                                 final boolean indexOnTheFly,
                                                 final boolean allowCompressed) {
        try {
            // we cannot merge compressed outputs, so don't compress if allowCompressed is false,
            // which is the case when we have a temporary output file for later merging
            if ( allowCompressed && stub.isCompressed() )
                stream = new BlockCompressedOutputStream(file);
            else
                stream = new PrintStream(new BufferedOutputStream(new FileOutputStream(file), BUFFER_SIZE));
        }
        catch(IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(file, "Unable to open target output stream", ex);
        }

        EnumSet<Options> options = stub.getWriterOptions(indexOnTheFly);
        VariantContextWriter writer = VariantContextWriterFactory.create(file, this.stream, stub.getMasterSequenceDictionary(), stub.getIndexCreator(), options);

        // if the stub says to test BCF, create a secondary writer to BCF and an 2 way out writer to send to both
        // TODO -- remove me when argument generateShadowBCF is removed
        if ( stub.alsoWriteBCFForTest() && ! VariantContextWriterFactory.isBCFOutput(file, options)) {
            final File bcfFile = BCF2Utils.shadowBCF(file);
            if ( bcfFile != null ) {
                FileOutputStream bcfStream;
                try {
                    bcfStream = new FileOutputStream(bcfFile);
                } catch (FileNotFoundException e) {
                    throw new RuntimeException(bcfFile + ": Unable to create BCF writer", e);
                }

                VariantContextWriter bcfWriter = VariantContextWriterFactory.create(bcfFile, bcfStream, stub.getMasterSequenceDictionary(), stub.getIndexCreator(), options);
                writer = new TestWriter(writer, bcfWriter);
            }
        }

        return writer;
    }

    /**
     * Check the return from PrintStream.checkError() if underlying stream is a java.io.PrintStream
     * @return true if PrintStream.checkError() returned true, false otherwise
     */
    public boolean checkError(){
        if ( stream instanceof PrintStream )
            return ((PrintStream) stream).checkError();
        return false;
    }

    private final static class TestWriter implements VariantContextWriter {
        final List<VariantContextWriter> writers;

        private TestWriter(final VariantContextWriter ... writers) {
            this.writers = Arrays.asList(writers);
        }

        @Override
        public void writeHeader(final VCFHeader header) {
            for ( final VariantContextWriter writer : writers ) writer.writeHeader(header);
        }

        @Override
        public void close() {
            for ( final VariantContextWriter writer : writers ) writer.close();
        }

        @Override
        public void add(final VariantContext vc) {
            for ( final VariantContextWriter writer : writers ) writer.add(vc);
        }

        /**
         * Check the return from PrintStream.checkError() if underlying stream for a java.io.PrintStream
         * @return false, no error since the underlying stream is not a java.io.PrintStream
         */
        public boolean checkError(){
            return false;
        }
    }

    public void add(VariantContext vc) {
        if ( closed ) throw new ReviewedGATKException("Attempting to write to a closed VariantContextWriterStorage " + vc.getStart() + " storage=" + this);
        writer.add(vc);
    }

    /**
     * initialize this VCF header
     *
     * @param header  the header
     */
    public void writeHeader(VCFHeader header) {
        writer.writeHeader(header);
    }

    /**
     * Close the VCF storage object.
     */
    public void close() {
        writer.close();
        closed = true;
    }

    public void mergeInto(VariantContextWriterStorage target) {
        try {
            if ( ! closed )
                throw new ReviewedGATKException("Writer not closed, but we are merging into the file!");
            final String targetFilePath = target.file != null ? target.file.getAbsolutePath() : "/dev/stdin";
            logger.debug(String.format("Merging VariantContextWriterStorage from %s into %s", file.getAbsolutePath(), targetFilePath));

            // use the feature manager to determine the right codec for the tmp file
            // that way we don't assume it's a specific type
            final FeatureManager.FeatureDescriptor fd = new FeatureManager().getByFiletype(file);
            if ( fd == null )
                throw new UserException.LocalParallelizationProblem(file);

            final FeatureCodec codec = fd.getCodec();
            final AbstractFeatureReader<Feature, ?> source = AbstractFeatureReader.getFeatureReader(file.getAbsolutePath(), codec, false);

            for ( final Feature vc : source.iterator() ) {
                target.writer.add((VariantContext) vc);
            }

            source.close();
            file.delete(); // this should be last to aid in debugging when the process fails
        } catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(file, "Error reading file in VCFWriterStorage: ", e);
        }
    }

}
