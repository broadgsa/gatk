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

import htsjdk.samtools.SAMFileWriter;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.utils.sam.GATKSAMFileWriter;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.io.OutputStream;
import java.lang.reflect.Type;

/**
 * Insert a SAMFileWriterStub  instead of a full-fledged concrete OutputStream implementations.
 */
public class SAMFileWriterArgumentTypeDescriptor extends ArgumentTypeDescriptor {

    /**
     * The engine into which output stubs should be fed.
     */
    private final GenomeAnalysisEngine engine;

    /**
     * The default location to which data should be written if the user specifies no such location.
     */
    private final OutputStream defaultOutputStream;

    /**
     * Create a new SAMFileWriter argument, notifying the given engine when that argument has been created.
     * @param engine Engine to add SAMFileWriter output to.
     * @param defaultOutputStream the target for the data
     */
    public SAMFileWriterArgumentTypeDescriptor( GenomeAnalysisEngine engine, OutputStream defaultOutputStream ) {
        this.engine = engine;
        this.defaultOutputStream = defaultOutputStream;
    }

    @Override
    public boolean supports( Class type ) {
        return SAMFileWriter.class.equals(type) || GATKSAMFileWriter.class.equals(type);
    }

    @Override
    public boolean createsTypeDefault(ArgumentSource source) {
        return !source.isRequired() && source.defaultsToStdout();
    }

    @Override
    public String typeDefaultDocString(ArgumentSource source) {
        return "stdout";
    }

    @Override
    public Object createTypeDefault(ParsingEngine parsingEngine,ArgumentSource source, Type type) {
        if(source.isRequired() || !source.defaultsToStdout())
            throw new ReviewedGATKException("BUG: tried to create type default for argument type descriptor that can't support a type default.");
        SAMFileWriterStub stub = new SAMFileWriterStub(engine,defaultOutputStream);
        engine.addOutput(stub);
        return stub;
    }

    @Override
    public Object parse( ParsingEngine parsingEngine, ArgumentSource source, Type type, ArgumentMatches matches )  {
        // Extract all possible parameters that could be passed to a BAM file writer?
        ArgumentDefinition bamArgumentDefinition = createDefaultArgumentDefinition(source);
        ArgumentMatchValue writerFileName = getArgumentValue( bamArgumentDefinition, matches );

        // Create the stub
        SAMFileWriterStub stub = null;      // stub = new SAMFileWriterStub(engine, defaultOutputStream);

        if (writerFileName != null &&  writerFileName.asFile() != null ) {
            stub = new SAMFileWriterStub(engine, writerFileName.asFile());

            // WARNING: Side effects required by engine!
            parsingEngine.addTags(stub,getArgumentTags(matches));
            engine.addOutput(stub);
        }

        return stub;
    }

}
