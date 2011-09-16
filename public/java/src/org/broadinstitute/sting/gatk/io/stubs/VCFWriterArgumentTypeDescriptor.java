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

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;
import java.io.OutputStream;
import java.lang.reflect.Type;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;

/**
 * Injects new command-line arguments into the system providing support for the genotype writer.
 *
 * @author mhanna
 * @version 0.1
 */
public class VCFWriterArgumentTypeDescriptor extends ArgumentTypeDescriptor {
    public static final String NO_HEADER_ARG_NAME = "NO_HEADER";
    public static final String SITES_ONLY_ARG_NAME = "sites_only";
    public static final HashSet<String> SUPPORTED_ZIPPED_SUFFIXES = new HashSet<String>();

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
        return VCFWriter.class.equals(type);
    }

    @Override
    public List<ArgumentDefinition> createArgumentDefinitions( ArgumentSource source ) {
        return Arrays.asList( createDefaultArgumentDefinition(source),createNoHeaderArgumentDefinition(),createSitesOnlyArgumentDefinition());
    }

    /**
     * This command-line argument descriptor does want to override the provided default value.
     * @return true always.
     */
    @Override
    public boolean createsTypeDefault(ArgumentSource source) {
        return source.isRequired();
    }

    @Override
    public String typeDefaultDocString(ArgumentSource source) {
        return "stdout";
    }

    @Override
    public Object createTypeDefault(ParsingEngine parsingEngine,ArgumentSource source, Type type) {
        if(!source.isRequired())
            throw new ReviewedStingException("BUG: tried to create type default for argument type descriptor that can't support a type default.");        
        VCFWriterStub stub = new VCFWriterStub(engine, defaultOutputStream, false, argumentSources, false, false);
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
        String writerFileName = getArgumentValue(defaultArgumentDefinition,matches);
        File writerFile = writerFileName != null ? new File(writerFileName) : null;

        // This parser has been passed a null filename and the GATK is not responsible for creating a type default for the object;
        // therefore, the user must have failed to specify a type default
        if(writerFile == null && !source.isRequired())
            throw new MissingArgumentValueException(defaultArgumentDefinition);

        // Should we compress the output stream?
        boolean compress = isCompressed(writerFileName);

        boolean skipWritingHeader = argumentIsPresent(createNoHeaderArgumentDefinition(),matches);
        boolean doNotWriteGenotypes = argumentIsPresent(createSitesOnlyArgumentDefinition(),matches);

        // Create a stub for the given object.
        VCFWriterStub stub = (writerFile != null) ? new VCFWriterStub(engine, writerFile, compress, argumentSources, skipWritingHeader, doNotWriteGenotypes)
                                                  : new VCFWriterStub(engine, defaultOutputStream, compress, argumentSources, skipWritingHeader, doNotWriteGenotypes);

        // WARNING: Side effects required by engine!
        parsingEngine.addTags(stub,getArgumentTags(matches));
        engine.addOutput(stub);

        return stub;
    }

    /**
     * Creates the optional compression level argument for the BAM file.
     * @return Argument definition for the BAM file itself.  Will not be null.
     */
    private ArgumentDefinition createNoHeaderArgumentDefinition() {
        return new ArgumentDefinition( ArgumentIOType.ARGUMENT,
                                       boolean.class,
                                       NO_HEADER_ARG_NAME,
                                       NO_HEADER_ARG_NAME,
                                       "Don't output the usual VCF header tag with the command line. FOR DEBUGGING PURPOSES ONLY. This option is required in order to pass integration tests.",
                                       false,
                                       true,
                                       false,
                                       true,
                                       null,
                                       null,
                                       null,
                                       null );
    }

    /**
     * Creates the optional compression level argument for the BAM file.
     * @return Argument definition for the BAM file itself.  Will not be null.
     */
    private ArgumentDefinition createSitesOnlyArgumentDefinition() {
        return new ArgumentDefinition( ArgumentIOType.ARGUMENT,
                                       boolean.class,
                                       SITES_ONLY_ARG_NAME,
                                       SITES_ONLY_ARG_NAME,
                                       "Just output sites without genotypes (i.e. only the first 8 columns of the VCF)",
                                       false,
                                       true,
                                       false,
                                       true,
                                       null,
                                       null,
                                       null,
                                       null );
    }

    /**
     * Returns true if the file will be compressed.
     * @param writerFileName Name of the file
     * @return true if the file will be compressed.
     */
    public static boolean isCompressed(String writerFileName) {
        return writerFileName != null && SUPPORTED_ZIPPED_SUFFIXES.contains(getFileSuffix(writerFileName));
    }

    /**
     * Returns a lower-cased version of the suffix of the provided file.
     * @param fileName the file name.  Must not be null.
     * @return lower-cased version of the file suffix.  Will not be null.
     */
    private static String getFileSuffix(String fileName) {
        int indexOfLastDot = fileName.lastIndexOf(".");
        if ( indexOfLastDot == -1 )
            return "";
        return fileName.substring(indexOfLastDot).toLowerCase();
    }

}
