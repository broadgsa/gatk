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

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.ProgressLoggerInterface;
import htsjdk.samtools.util.RuntimeIOException;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.io.stubs.SAMFileWriterStub;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.sam.SimplifyingSAMFileWriter;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

/**
 * Provides temporary storage for SAMFileWriters.
 *
 * @author mhanna
 * @version 0.1
 */
public class SAMFileWriterStorage implements SAMFileWriter, Storage<SAMFileWriter> {
    private final File file;
    private File referenceFasta;
    private SAMFileWriter writer;

    private static Logger logger = Logger.getLogger(SAMFileWriterStorage.class);

    public SAMFileWriterStorage( SAMFileWriterStub stub ) {
        this(stub,stub.getOutputFile());
    }

    public SAMFileWriterStorage( SAMFileWriterStub stub, File file ) {
        this.referenceFasta = stub.getReferenceFile();
        this.file = file;
        SAMFileWriterFactory factory = new SAMFileWriterFactory();
        // Enable automatic index creation for pre-sorted BAMs.
        if (stub.getFileHeader().getSortOrder().equals(SAMFileHeader.SortOrder.coordinate) && stub.getIndexOnTheFly())
            factory.setCreateIndex(true);
        if (stub.getGenerateMD5())
            factory.setCreateMd5File(true);
        // Adjust max records in RAM.
        // TODO -- this doesn't actually work because of a bug in Picard; do not use until fixed
        if(stub.getMaxRecordsInRam() != null)
            factory.setMaxRecordsInRam(stub.getMaxRecordsInRam());

        if(stub.getOutputFile() != null) {
            try {
                if (stub.getOutputFile().getName().toLowerCase().endsWith(".cram")) {
                    this.writer = createCRAMWriter(factory, stub.getFileHeader(), file, this.referenceFasta);
                } else {
                    this.writer = createBAMWriter(factory,stub.getFileHeader(),stub.isPresorted(),file,stub.getCompressionLevel());
                }
            } catch(RuntimeIOException ex) {
                throw new UserException.CouldNotCreateOutputFile(file,"file could not be created",ex);
            }
        }
        else if(stub.getOutputStream() != null){
            this.writer = factory.makeSAMWriter( stub.getFileHeader(), stub.isPresorted(), stub.getOutputStream());
        }
        else
            throw new UserException("Unable to write to SAM file; neither a target file nor a stream has been specified");

        // if we want to send the BAM file through the simplifying writer, wrap it here
        if ( stub.simplifyBAM() ) {
            this.writer = new SimplifyingSAMFileWriter(this.writer);
        }
    }

    public SAMFileHeader getFileHeader() {
        return writer.getFileHeader();
    }

    public void addAlignment( SAMRecord read ) {
        writer.addAlignment(read);
    }

    public void close() {
        try {
            writer.close();
        } catch (RuntimeIOException e) {
            throw new UserException.ErrorWritingBamFile(e.getMessage());
        }
    }

    public void mergeInto( SAMFileWriter targetStream ) {
        SAMFileReader reader = new SAMFileReader( file );
        try {
            CloseableIterator<SAMRecord> iterator = reader.iterator();
            while( iterator.hasNext() )
                targetStream.addAlignment( iterator.next() );
            iterator.close();
        }
        finally {
            reader.close();
            file.delete();
        }
    }

    private SAMFileWriter createCRAMWriter(final SAMFileWriterFactory factory,
                                           final SAMFileHeader header,
                                           final File file,
                                           final File referenceFasta) {
        return factory.makeCRAMWriter(header, file, referenceFasta);
    }

    private SAMFileWriter createBAMWriter(final SAMFileWriterFactory factory,
                                 final SAMFileHeader header,
                                 final boolean presorted,
                                 final File outputFile,
                                 final Integer compressionLevel) {
        SAMFileWriter writer;
        if(compressionLevel != null)
            writer = factory.makeBAMWriter(header, presorted, outputFile, compressionLevel);
        else
            writer = factory.makeBAMWriter(header, presorted, outputFile);

        // mhanna - 1 Mar 2011 - temporary hack until Picard generates an index file for empty BAMs --
        //                     - do a pre-initialization of the BAM file.
        try {
            Method prepareToWriteAlignmentsMethod = writer.getClass().getDeclaredMethod("prepareToWriteAlignments");
            if(prepareToWriteAlignmentsMethod != null) {
                prepareToWriteAlignmentsMethod.setAccessible(true);
                prepareToWriteAlignmentsMethod.invoke(writer);
            }
        }
        catch(NoSuchMethodException ex) {
            logger.info("Unable to call prepareToWriteAlignments method; this should be reviewed when Picard is updated.");
        }
        catch(IllegalAccessException ex) {
            logger.info("Unable to access prepareToWriteAlignments method; this should be reviewed when Picard is updated.");
        }
        catch(InvocationTargetException ex) {
            logger.info("Unable to invoke prepareToWriteAlignments method; this should be reviewed when Picard is updated.");
        }

        return writer;
    }

    @Override
    public void setProgressLogger(final ProgressLoggerInterface logger) {
        writer.setProgressLogger(logger);
    }
}
