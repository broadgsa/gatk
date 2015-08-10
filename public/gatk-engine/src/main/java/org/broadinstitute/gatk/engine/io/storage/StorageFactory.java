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
import org.broadinstitute.gatk.engine.io.stubs.SAMFileWriterStub;
import org.broadinstitute.gatk.engine.io.stubs.Stub;
import org.broadinstitute.gatk.engine.io.stubs.VariantContextWriterStub;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.io.File;

/**
 * Construct storage of the required type.
 *
 * @author mhanna
 * @version 0.1
 */
public class StorageFactory {
    /**
     * Disable storage factory construction.
     */
    private StorageFactory() {}

    /**
     * Gets the output storage associated with a given stub.
     * @param stub The stub for which to find / create the right output stream.
     * @param <T> Type of the stream to create.
     * @return Storage object with a facade of type T.
     */
    public static <T> Storage<T> createStorage( Stub<T> stub ) {
        return createStorage( stub, null );
    }

    /**
     * Gets the output storage associated with a given stub.
     * @param stub The stub for which to find / create the right output stream.
     * @param file The filename to which to write the file.
     * @param <T> Type of the stream to create.
     * @return Storage object with a facade of type T.
     */
     public static <T> Storage<T> createStorage( Stub<T> stub, File file ) {
        Storage storage;

        if(stub instanceof OutputStreamStub) {
            if( file != null )
                storage = new OutputStreamStorage((OutputStreamStub)stub,file);
            else
                storage = new OutputStreamStorage((OutputStreamStub)stub);
        }
        else if(stub instanceof SAMFileWriterStub) {
            if( file != null )
                storage = new SAMFileWriterStorage((SAMFileWriterStub)stub,file);
            else
                storage = new SAMFileWriterStorage((SAMFileWriterStub)stub);
        }
        else if(stub instanceof VariantContextWriterStub) {
            VariantContextWriterStub vcfWriterStub = (VariantContextWriterStub)stub;
            if( file != null )
                storage = new VariantContextWriterStorage(vcfWriterStub,file);
            else
                storage = new VariantContextWriterStorage(vcfWriterStub);
        }
        else
            throw new ReviewedGATKException("Unsupported stub type: " + stub.getClass().getName());

        return storage;
    }
}
