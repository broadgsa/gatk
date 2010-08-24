/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.gatk.io.stubs;

import org.broad.tribble.vcf.VCFWriter;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;

import java.io.File;
import java.io.OutputStream;
import java.util.List;
import java.util.Arrays;
import java.util.HashSet;

/**
 * Injects new command-line arguments into the system providing support for the genotype writer.
 *
 * @author mhanna
 * @version 0.1
 */
public class VCFWriterArgumentTypeDescriptor extends ArgumentTypeDescriptor {

    private static final HashSet<String> SUPPORTED_ZIPPED_SUFFIXES = new HashSet<String>();
    //
    // static list of zipped suffixes supported by this system.
    //
    static {
        SUPPORTED_ZIPPED_SUFFIXES.add(".gz");
        SUPPORTED_ZIPPED_SUFFIXES.add(".gzip");
    }

    /**
     * The engine into which output stubs should be fed.
     */
    private GenomeAnalysisEngine engine;

    /**
      * The default location to which data should be written if the user specifies no such location.
      */
    private final OutputStream defaultOutputStream;

    /**
     * Create a new GenotypeWriter argument, notifying the given engine when that argument has been created.
     * @param engine the engine to be notified.
     * @param defaultOutputStream the default output stream to be written to if nothing else is specified.
     */
    public VCFWriterArgumentTypeDescriptor(GenomeAnalysisEngine engine, OutputStream defaultOutputStream) {
        this.engine = engine;
        this.defaultOutputStream = defaultOutputStream;
    }

    /**
     * Reports whether this ArgumentTypeDescriptor supports the given type.
     * @param type The type to check.
     * @return True if the argument is a GenotypeWriter.
     */
    @Override
    public boolean supports( Class type ) {
        return VCFWriter.class.equals(type);
    }

    @Override
    public List<ArgumentDefinition> createArgumentDefinitions( ArgumentSource source ) {
        return Arrays.asList( createDefaultArgumentDefinition(source) );
    }

    /**
     * This command-line argument descriptor does want to override the provided default value.
     * @return true always.
     */
    @Override
    public boolean createsTypeDefault(ArgumentSource source,Class type) {
        return true;
    }

    @Override
    public Object createTypeDefault(ArgumentSource source,Class type) {
        VCFWriterStub stub = new VCFWriterStub(engine, defaultOutputStream, false);
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
    public Object parse( ArgumentSource source, Class type, ArgumentMatches matches )  {
        // Get the filename for the genotype file, if it exists.  If not, we'll need to send output to out.
        String writerFileName = getArgumentValue(createDefaultArgumentDefinition(source),matches);
        File writerFile = writerFileName != null ? new File(writerFileName) : null;

        // Should we compress the output stream?
        boolean compress = writerFileName != null && SUPPORTED_ZIPPED_SUFFIXES.contains(getFileSuffix(writerFileName));

        // Create a stub for the given object.
        VCFWriterStub stub = (writerFile != null) ? new VCFWriterStub(engine, writerFile, compress) : new VCFWriterStub(engine, System.out, compress);

        engine.addOutput(stub);

        return stub;
    }

    /**
     * Returns a lower-cased version of the suffix of the provided file.
     * @param fileName the file name.  Must not be null.
     * @return lower-cased version of the file suffix.  Will not be null.
     */
    private String getFileSuffix(String fileName) {
        int indexOfLastDot = fileName.lastIndexOf(".");
        if ( indexOfLastDot == -1 )
            return "";
        return fileName.substring(indexOfLastDot).toLowerCase();
    }
}
