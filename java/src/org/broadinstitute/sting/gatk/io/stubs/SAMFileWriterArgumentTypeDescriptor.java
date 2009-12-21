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

package org.broadinstitute.sting.gatk.io.stubs;

import org.broadinstitute.sting.utils.cmdLine.*;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import net.sf.samtools.SAMFileWriter;

import java.util.List;
import java.util.Arrays;
import java.io.File;

/**
 * Insert a SAMFileWriterStub  instead of a full-fledged concrete OutputStream implementations.
 */
public class SAMFileWriterArgumentTypeDescriptor extends ArgumentTypeDescriptor {
    private static final String DEFAULT_ARGUMENT_FULLNAME = "outputBAM";
    private static final String DEFAULT_ARGUMENT_SHORTNAME = "ob";

    private static final String COMPRESSION_FULLNAME = "bam_compression";
    private static final String COMPRESSION_SHORTNAME = "compress";

    /**
     * The engine into which output stubs should be fed.
     */
    private GenomeAnalysisEngine engine;

    /**
     * Create a new SAMFileWriter argument, notifying the given engine when that argument has been created.
     * @param engine Engine to add SAMFileWriter output to.
     */
    public SAMFileWriterArgumentTypeDescriptor( GenomeAnalysisEngine engine ) {
        this.engine = engine;
    }
    

    @Override
    public boolean supports( Class type ) {
        return SAMFileWriter.class.isAssignableFrom(type);
    }

    @Override
    public List<ArgumentDefinition> createArgumentDefinitions( ArgumentSource source ) {
        return Arrays.asList( createBAMArgumentDefinition(source),
                              createBAMCompressionArgumentDefinition(source) );
    }

    @Override
    public Object parse( ArgumentSource source, Class type, ArgumentMatches matches )  {
        String writerFileName = getArgumentValue( createBAMArgumentDefinition(source), matches );
        if( writerFileName == null )
            throw new StingException("SAM file compression was supplied, but no associated writer was supplied with it.");

        SAMFileWriterStub stub = new SAMFileWriterStub(engine, new File(writerFileName));

        String compressionLevelText = getArgumentValue( createBAMCompressionArgumentDefinition(source), matches );
        Integer compressionLevel = compressionLevelText != null ? Integer.valueOf(compressionLevelText) : null;
        if( compressionLevel != null )
            stub.setCompressionLevel(compressionLevel);

        engine.addOutput(stub);

        return stub;
    }

    /**
     * Gets the definition of the argument representing the BAM file itself.
     * @param source Argument source for the BAM file.  Must not be null.
     * @return Argument definition for the BAM file itself.  Will not be null.
     */
    private ArgumentDefinition createBAMArgumentDefinition(ArgumentSource source) {
        Argument description = this.getArgumentDescription(source);

        boolean isFullNameProvided = description.fullName().trim().length() > 0;
        boolean isShortNameProvided = description.shortName().trim().length() > 0;

        String fullName = isFullNameProvided ? description.fullName().trim() : DEFAULT_ARGUMENT_FULLNAME;

        // If the short name is provided, use that.  If the user hasn't provided any names at all, use
        // the default.  If somewhere in the middle, leave the short name blank.
        String shortName;
        if( isShortNameProvided )
            shortName = description.shortName().trim();
        else if( !isFullNameProvided )
            shortName = DEFAULT_ARGUMENT_SHORTNAME;
        else
            shortName = null;

        return new ArgumentDefinition( fullName,
                                       shortName,
                                       getDoc(source),
                                       isRequired(source),
                                       false,
                                       source.isMultiValued(),
                                       getExclusiveOf(source),
                                       getValidationRegex(source) );
    }

    /**
     * Creates the optional compression level argument for the BAM file.
     * @param source Argument source for the BAM file.  Must not be null.
     * @return Argument definition for the BAM file itself.  Will not be null.
     */
    private ArgumentDefinition createBAMCompressionArgumentDefinition(ArgumentSource source) {
        return new ArgumentDefinition( COMPRESSION_FULLNAME,
                                       COMPRESSION_SHORTNAME,
                                       "Compression level to use for writing BAM files",
                                       false,
                                       false,
                                       false,
                                       null,
                                       null );
    }
}
