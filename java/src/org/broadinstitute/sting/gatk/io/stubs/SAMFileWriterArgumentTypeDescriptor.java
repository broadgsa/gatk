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
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import net.sf.samtools.SAMFileWriter;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.lang.annotation.Annotation;
import java.util.List;
import java.util.Arrays;
import java.io.File;
import java.io.OutputStream;

/**
 * Insert a SAMFileWriterStub  instead of a full-fledged concrete OutputStream implementations.
 */
public class SAMFileWriterArgumentTypeDescriptor extends ArgumentTypeDescriptor {
    private static final String DEFAULT_ARGUMENT_FULLNAME = "outputBAM";
    private static final String DEFAULT_ARGUMENT_SHORTNAME = "ob";

    private static final String COMPRESSION_FULLNAME = "bam_compression";
    private static final String COMPRESSION_SHORTNAME = "compress";

    private static final String CREATE_INDEX_FULLNAME = "index_output_bam_on_the_fly";

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
        return SAMFileWriter.class.equals(type) || StingSAMFileWriter.class.equals(type);
    }

    @Override
    public List<ArgumentDefinition> createArgumentDefinitions( ArgumentSource source ) {
        return Arrays.asList( createBAMArgumentDefinition(source),
                              createBAMCompressionArgumentDefinition(source),
                              createWriteIndexArgumentDefinition(source));
    }

    @Override
    public boolean createsTypeDefault(ArgumentSource source) {
        return source.isRequired();
    }

    @Override
    public Object createTypeDefault(ParsingEngine parsingEngine,ArgumentSource source) {
        if(!source.isRequired())
            throw new ReviewedStingException("BUG: tried to create type default for argument type descriptor that can't support a type default.");
        SAMFileWriterStub stub = new SAMFileWriterStub(engine,defaultOutputStream);
        engine.addOutput(stub);
        return stub;
    }

    @Override
    public Object parse( ParsingEngine parsingEngine, ArgumentSource source, Class type, ArgumentMatches matches )  {
        ArgumentDefinition bamArgumentDefinition = createBAMArgumentDefinition(source);
        String writerFileName = getArgumentValue( bamArgumentDefinition, matches );

        // This parser has been passed a null filename and the GATK is not responsible for creating a type default for the object;
        // therefore, the user must have failed to specify a type default
        if(writerFileName == null && !source.isRequired())
            throw new MissingArgumentValueException(bamArgumentDefinition);

        SAMFileWriterStub stub = new SAMFileWriterStub(engine, new File(writerFileName));

        String compressionLevelText = getArgumentValue( createBAMCompressionArgumentDefinition(source), matches );
        Integer compressionLevel = compressionLevelText != null ? Integer.valueOf(compressionLevelText) : null;
        if( compressionLevel != null )
            stub.setCompressionLevel(compressionLevel);

        stub.setIndexOnTheFly(argumentIsPresent(createWriteIndexArgumentDefinition(source),matches));

        // WARNING: Side effects required by engine!
        parsingEngine.addTags(stub,getArgumentTags(matches));
        engine.addOutput(stub);

        return stub;
    }

    /**
     * Gets the definition of the argument representing the BAM file itself.
     * @param source Argument source for the BAM file.  Must not be null.
     * @return Argument definition for the BAM file itself.  Will not be null.
     */
    private ArgumentDefinition createBAMArgumentDefinition(ArgumentSource source) {
        Annotation annotation = getArgumentAnnotation(source);
        return new ArgumentDefinition( annotation,
                                       ArgumentIOType.getIOType(annotation),
                                       source.field.getType(),
                                       DEFAULT_ARGUMENT_FULLNAME,
                                       DEFAULT_ARGUMENT_SHORTNAME,
                                       ArgumentDefinition.getDoc(annotation),
                                       false,
                                       false,
                                       source.isMultiValued(),
                                       source.isHidden(),
                                       null,
                                       null,
                                       null,
                                       null);
    }

    /**
     * Creates the optional compression level argument for the BAM file.
     * @param source Argument source for the BAM file.  Must not be null.
     * @return Argument definition for the BAM file itself.  Will not be null.
     */
    private ArgumentDefinition createBAMCompressionArgumentDefinition(ArgumentSource source) {
        return new ArgumentDefinition( ArgumentIOType.ARGUMENT,
                                       int.class,
                                       COMPRESSION_FULLNAME,
                                       COMPRESSION_SHORTNAME,
                                       "Compression level to use for writing BAM files",
                                       false,
                                       false,
                                       false,
                                       source.isHidden(),
                                       null,
                                       null,
                                       null,
                                       null );
    }

    private ArgumentDefinition createWriteIndexArgumentDefinition(ArgumentSource source) {
        return new ArgumentDefinition( ArgumentIOType.ARGUMENT,
                                       boolean.class,
                                       CREATE_INDEX_FULLNAME,
                                       null,
                                       "Create a BAM index on-the-fly while writing the resulting file.",
                                       false,
                                       false,
                                       false,
                                       source.isHidden(),
                                       null,
                                       null,
                                       null,
                                       null );
    }
}
