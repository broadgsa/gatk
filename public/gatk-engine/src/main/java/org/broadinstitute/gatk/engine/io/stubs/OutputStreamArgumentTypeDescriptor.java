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
import org.broadinstitute.gatk.utils.exceptions.DynamicClassResolutionException;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.io.File;
import java.io.OutputStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.Type;

/**
 * Insert an OutputStreamStub instead of a full-fledged concrete OutputStream implementations.
 */
public class OutputStreamArgumentTypeDescriptor extends ArgumentTypeDescriptor {
    /**
     * The engine into which output stubs should be fed.
     */
    private final GenomeAnalysisEngine engine;

    /**
     * The default output stream to write to write this info if
     */
    private final OutputStream defaultOutputStream;

    /**
     * Create a new OutputStream argument, notifying the given engine when that argument has been created.
     * @param engine Engine to add SAMFileWriter output to.
     * @param defaultOutputStream Default target for output file.
     */
    public OutputStreamArgumentTypeDescriptor(GenomeAnalysisEngine engine,OutputStream defaultOutputStream) {
        this.engine = engine;
        this.defaultOutputStream = defaultOutputStream;
    }

    @Override
    public boolean supports( Class type ) {
        return getConstructorForClass(type) != null;
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
        OutputStreamStub stub = new OutputStreamStub(defaultOutputStream);
        engine.addOutput(stub);
        return createInstanceOfClass((Class)type,stub);
    }

    @Override
    public Object parse( ParsingEngine parsingEngine, ArgumentSource source, Type type, ArgumentMatches matches )  {
        ArgumentDefinition definition = createDefaultArgumentDefinition(source);
        String fileName = getArgumentValue( definition, matches ).asString();

        // This parser has been passed a null filename and the GATK is not responsible for creating a type default for the object;
        // therefore, the user must have failed to specify a type default
        if(fileName == null && source.isRequired())
            throw new MissingArgumentValueException(definition);

        OutputStreamStub stub = new OutputStreamStub(new File(fileName));

        engine.addOutput(stub);

        Object result = createInstanceOfClass(makeRawTypeIfNecessary(type),stub);
        // WARNING: Side effects required by engine!
        parsingEngine.addTags(result,getArgumentTags(matches));
        
        return result;
    }

    /**
     * Retrieves the constructor for an object that takes exactly one argument: an output stream.
     * @param type Type for which to go constructor spelunking.
     * @return Constructor, if available.  Null, if not.
     */
    private Constructor<OutputStream> getConstructorForClass( Class type ) {
        try {
            return type.getConstructor( OutputStream.class );
        }
        catch( NoSuchMethodException ex ) {
            return null;
        }
    }

    /**
     * Creat a new instance of the class accepting a single outputstream constructor.
     * @param type Type of object to create.
     * @param outputStream resulting output stream.
     * @return A new instance of the outputstream-derived class.
     */
    private Object createInstanceOfClass(Class type,OutputStream outputStream) {
        try {
            return getConstructorForClass(type).newInstance(outputStream);
        } catch (Exception e) {
            throw new DynamicClassResolutionException(type, e);
        }
    }
}
