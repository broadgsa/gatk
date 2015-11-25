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

import org.broadinstitute.gatk.engine.io.stubs.OutputStreamStub;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.io.*;
import java.nio.channels.Channels;
import java.nio.channels.FileChannel;
import java.nio.channels.WritableByteChannel;

public class OutputStreamStorage extends OutputStream implements Storage<OutputStream> {
    /**
     * File to which data will temporarily be written.
     */
    private final File file;

    /**
     * Stream to which data in this shard will be written.
     */
    private final OutputStream outputStream;

    /**
     * Create a new storage area with the given stub.
     * @param stub
     */
    public OutputStreamStorage( OutputStreamStub stub ) {
        if( stub.getOutputFile() != null ) {
            this.file = stub.getOutputFile();
            this.outputStream = initializeOutputStream(stub.getOutputFile());
        }
        else if( stub.getOutputStream() != null ) {
            this.file = null;
            this.outputStream = stub.getOutputStream();           
        }
        else
            throw new ReviewedGATKException("Not enough information to create storage for an OutputStream; need either a file or an existing output stream");
    }

    public OutputStreamStorage( OutputStreamStub stub, File file ) {
        this.file = file;
        this.outputStream = initializeOutputStream(file);
    }

    private OutputStream initializeOutputStream( File file ) {
        try {
            return new FileOutputStream( file );
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotCreateOutputFile(file, "Unable to open output stream for file", ex);
        }
    }

    /**
     * @{inheritDoc}
     */
    public void flush() throws IOException {
        outputStream.flush();
    }

    /**
     * @{inheritDoc}
     */
    public void close() {
        // Don't close System.out or System.err; this'll cause trouble
        // with subsequent code running in this VM.
        if( outputStream == System.out || outputStream == System.err )
            return;
        
        try {
            outputStream.close();
        }
        catch( IOException ex ) {
            throw new UserException.CouldNotCreateOutputFile(file, "Unable to close output stream", ex );
        }
    }

    /**
     * @{inheritDoc}
     */
    public void write( byte[] b ) throws IOException {
        outputStream.write(b);
    }

    /**
     * @{inheritDoc}
     */
    public void write( byte[] b, int off, int len ) throws IOException {
        outputStream.write(b, off, len);
    }

    /**
     * @{inheritDoc}
     */
    public void write( int b ) throws IOException {
        outputStream.write(b);
    }


    public void mergeInto( OutputStream targetStream ) {
        FileInputStream sourceStream = null;
        try {
            sourceStream = new FileInputStream( file );
            FileChannel sourceChannel = sourceStream.getChannel();

            WritableByteChannel targetChannel = Channels.newChannel( targetStream );
            sourceChannel.transferTo( 0, sourceChannel.size(), targetChannel );

            sourceStream.close();
            file.delete();
        }
        catch( FileNotFoundException ex ) {
            throw new UserException.CouldNotReadInputFile(file, "Unable to open input stream for file", ex);
        }
        catch( IOException ex ) {
            throw new UserException.CouldNotReadInputFile(file, "Unable to transfer contents of file", ex);
        }
    }
}
