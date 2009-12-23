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

package org.broadinstitute.sting.gatk.io.storage;

import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;

import java.io.*;

import org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterStub;

/**
 * Provides temporary storage for SAMFileWriters.
 *
 * @author mhanna
 * @version 0.1
 */
public class SAMFileWriterStorage implements SAMFileWriter, Storage<SAMFileWriter> {
    private final File file;
    private final SAMFileWriter writer;

    public SAMFileWriterStorage( SAMFileWriterStub stub ) {
        this(stub,stub.getSAMFile());   
    }

    public SAMFileWriterStorage( SAMFileWriterStub stub, File file ) {
        this.file = file;
        if( stub.getCompressionLevel() != null )
            this.writer = new SAMFileWriterFactory().makeBAMWriter( stub.getSAMFileHeader(), true, file, stub.getCompressionLevel() );
        else
            this.writer = new SAMFileWriterFactory().makeBAMWriter( stub.getSAMFileHeader(), true, file );
    }

    public void addAlignment( SAMRecord read ) {
        writer.addAlignment(read);
    }

    public SAMFileHeader getFileHeader() {
        return writer.getFileHeader();
    }

    public void close() {
        writer.close();
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

}