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

import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.io.File;
import java.io.OutputStream;
import java.lang.reflect.Type;
import java.util.Collection;

/**
 * Injects new command-line arguments into the system providing support for the genotype writer.
 *
 * @author mhanna
 * @version 0.1
 */
public class VCFWriterArgumentTypeDescriptor extends ArgumentTypeDescriptor {

    /**
     * The engine into which output stubs should be fed.
     */
    private final GenomeAnalysisEngine engine;

    /**
      * The default location to which data should be written if the user specifies no such location.
      */
    private final OutputStream defaultOutputStream;

    /**
     * The sources into which arguments were injected.
     */
    private final Collection<Object> argumentSources;

    /**
     * Create a new GenotypeWriter argument, notifying the given engine when that argument has been created.
     * @param engine the engine to be notified.
     * @param defaultOutputStream the default output stream to be written to if nothing else is specified.
     * @param argumentSources sources from which command-line arguments should be derived.
     */
    public VCFWriterArgumentTypeDescriptor(GenomeAnalysisEngine engine, OutputStream defaultOutputStream, Collection<Object> argumentSources) {
        this.engine = engine;
        this.defaultOutputStream = defaultOutputStream;
        this.argumentSources = argumentSources;
    }

    /**
     * Reports whether this ArgumentTypeDescriptor supports the given type.
     * @param type The type to check.
     * @return True if the argument is a GenotypeWriter.
     */
    @Override
    public boolean supports( Class type ) {
        return VariantContextWriter.class.equals(type);
    }

    /**
     * This command-line argument descriptor does want to override the provided default value.
     * @return true always.
     */
    @Override
    public boolean createsTypeDefault(ArgumentSource source) {
        return !source.isRequired() && source.defaultsToStdout();
    }

    @Override
    public String typeDefaultDocString(ArgumentSource source) {
        return "stdout";
    }

    @Override
    public Object createTypeDefault(ParsingEngine parsingEngine, ArgumentSource source, Type type) {
        if(source.isRequired() || !source.defaultsToStdout())
            throw new ReviewedGATKException("BUG: tried to create type default for argument type descriptor that can't support a type default.");        
        VariantContextWriterStub stub = new VariantContextWriterStub(engine, defaultOutputStream, argumentSources);
        engine.addOutput(stub);
        return stub;
    }

    /**
     * Convert the given argument matches into a single object suitable for feeding into the ArgumentSource.
     * @param source Source for this argument.
     * @param type not used
     * @param matches Matches that match with this argument.
     * @return Transform from the matches into the associated argument.
     */
    @Override
    public Object parse( ParsingEngine parsingEngine, ArgumentSource source, Type type, ArgumentMatches matches )  {
        ArgumentDefinition defaultArgumentDefinition = createDefaultArgumentDefinition(source);
        // Get the filename for the genotype file, if it exists.  If not, we'll need to send output to out.
        ArgumentMatchValue writerFileName = getArgumentValue(defaultArgumentDefinition,matches);
        File writerFile = writerFileName != null ? writerFileName.asFile() : null;

        // This parser has been passed a null filename and the GATK is not responsible for creating a type default for the object;
        // therefore, the user must have failed to specify a type default
        if(writerFile == null && source.isRequired())
            throw new MissingArgumentValueException(defaultArgumentDefinition);

        // Create a stub for the given object.
        final VariantContextWriterStub stub = (writerFile != null)
                ? new VariantContextWriterStub(engine, writerFile, argumentSources)
                : new VariantContextWriterStub(engine, defaultOutputStream, argumentSources);

        stub.setCompressed(isCompressed(writerFileName == null ? null: writerFileName.asString()));

        // WARNING: Side effects required by engine!
        parsingEngine.addTags(stub,getArgumentTags(matches));
        engine.addOutput(stub);

        return stub;
    }
}
