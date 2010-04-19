/*
 * Copyright (c) 2010 The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the ”Software”), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED ”AS IS”, WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.io.stubs;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;

import java.io.OutputStream;
import java.io.File;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;

/**
 * Insert an OutputStreamStub instead of a full-fledged concrete OutputStream implementations.
 */
public class OutputStreamArgumentTypeDescriptor extends ArgumentTypeDescriptor {
    /**
     * The engine into which output stubs should be fed.
     */
    private GenomeAnalysisEngine engine;

    /**
     * Create a new OutputStream argument, notifying the given engine when that argument has been created.
     * @param engine Engine to add SAMFileWriter output to.
     */
    public OutputStreamArgumentTypeDescriptor( GenomeAnalysisEngine engine ) {
        this.engine = engine;    
    }

    @Override
    public boolean supports( Class type ) {
        return getConstructorForClass(type) != null;
    }

    @Override
    public Object parse( ArgumentSource source, Class type, ArgumentMatches matches )  {
        ArgumentDefinition definition = createDefaultArgumentDefinition(source);
        String fileName = getArgumentValue( definition, matches );

        OutputStreamStub stub = new OutputStreamStub(new File(fileName));

        engine.addOutput(stub);

        try {
            return getConstructorForClass(type).newInstance(stub);
        }
        catch( InstantiationException ex ) {
            throw new StingException("Could not instantiate class with OutputStream constructor: " + type.getName());
        }
        catch( IllegalAccessException ex ) {
            throw new StingException("Could not access class with OutputStream constructor: " + type.getName());                        
        }
        catch( InvocationTargetException ex ) {
            throw new StingException("Could not invoke constructor for class with OutputStream constructor: " + type.getName());                        
        }
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
}
