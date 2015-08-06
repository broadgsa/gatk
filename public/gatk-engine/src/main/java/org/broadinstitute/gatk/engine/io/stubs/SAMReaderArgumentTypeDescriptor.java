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

import htsjdk.samtools.SAMFileReader;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.sam.SAMReaderBuilder;

import java.lang.reflect.Type;

/**
 * Describe how to parse SAMReaders.
 */
public class SAMReaderArgumentTypeDescriptor extends ArgumentTypeDescriptor {
    /**
     * The engine into which output stubs should be fed.
     */
    private GenomeAnalysisEngine engine;

    /**
     * Create a new SAMFileReader argument, notifying the given engine when that argument has been created.
     * @param engine engine
     */
    public SAMReaderArgumentTypeDescriptor(GenomeAnalysisEngine engine) {
        this.engine = engine;
    }

    @Override
    public boolean supports( Class type ) {
        return SAMFileReader.class.isAssignableFrom(type);
    }

    @Override
    public Object parse( ParsingEngine parsingEngine, ArgumentSource source, Type type, ArgumentMatches matches ) {
        SAMReaderBuilder builder = new SAMReaderBuilder();

        ArgumentMatchValue readerFileName = getArgumentValue( createDefaultArgumentDefinition(source), matches );

        if( readerFileName == null )
            throw new UserException.CommandLineException("SAM file compression was supplied, but no associated writer was supplied with it.");

        builder.setSAMFile(readerFileName.asFile());

        // WARNING: Skipping required side-effect because stub is impossible to generate.
        engine.addInput(source, builder);

        // MASSIVE KLUDGE!  SAMFileReader is tricky to implement and we don't yet have a stub.  Return null, then
        // let the output tracker load it in.
        // TODO: Add a stub for SAMReader.
        return null;
    }
}
