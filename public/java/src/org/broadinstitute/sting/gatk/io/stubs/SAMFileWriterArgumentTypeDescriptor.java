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

import net.sf.samtools.SAMFileWriter;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;
import java.io.OutputStream;
import java.lang.annotation.Annotation;
import java.lang.reflect.Type;
import java.util.Arrays;
import java.util.List;

/**
 * Insert a SAMFileWriterStub  instead of a full-fledged concrete OutputStream implementations.
 */
public class SAMFileWriterArgumentTypeDescriptor extends ArgumentTypeDescriptor {
    public static final String DEFAULT_ARGUMENT_FULLNAME = "outputBAM";
    public static final String DEFAULT_ARGUMENT_SHORTNAME = "ob";

    public static final String COMPRESSION_FULLNAME = "bam_compression";
    public static final String COMPRESSION_SHORTNAME = "compress";

    public static final String SIMPLIFY_BAM_FULLNAME = "simplifyBAM";
    public static final String SIMPLIFY_BAM_SHORTNAME = SIMPLIFY_BAM_FULLNAME;

    public static final String DISABLE_INDEXING_FULLNAME = "disable_bam_indexing";
    public static final String ENABLE_MD5_FULLNAME = "generate_md5";

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
                              disableWriteIndexArgumentDefinition(source),
                              enableMD5GenerationArgumentDefinition(source),
                              createSimplifyBAMArgumentDefinition(source));
    }

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
        SAMFileWriterStub stub = new SAMFileWriterStub(engine,defaultOutputStream);
        engine.addOutput(stub);
        return stub;
    }

    @Override
    public Object parse( ParsingEngine parsingEngine, ArgumentSource source, Type type, ArgumentMatches matches )  {
        // Extract all possible parameters that could be passed to a BAM file writer?
        ArgumentDefinition bamArgumentDefinition = createBAMArgumentDefinition(source);
        String writerFileName = getArgumentValue( bamArgumentDefinition, matches );

        String compressionLevelText = getArgumentValue( createBAMCompressionArgumentDefinition(source), matches );
        Integer compressionLevel = compressionLevelText != null ? Integer.valueOf(compressionLevelText) : null;

        Boolean indexOnTheFly = !argumentIsPresent(disableWriteIndexArgumentDefinition(source),matches) ? true : null;
        Boolean generateMD5 = argumentIsPresent(this.enableMD5GenerationArgumentDefinition(source),matches) ? true : null;
        Boolean simplifyBAM = argumentIsPresent(createSimplifyBAMArgumentDefinition(source),matches);

        // Validate the combination of parameters passed in.

        // This parser has been passed a null filename and the GATK is not responsible for creating a type default for the object;
        // therefore, the user must have failed to specify a type default
        if(writerFileName == null) {
            if(!source.isRequired())
                throw new MissingArgumentValueException(bamArgumentDefinition);
            if(generateMD5)
                throw new ArgumentException("MD5 generation specified, but no output file specified.  If md5 generation is desired, please specify a BAM output file and an md5 file will be written alongside.");
        }

        // Create the stub and set parameters.
        SAMFileWriterStub stub = new SAMFileWriterStub(engine, new File(writerFileName));

        if( compressionLevel != null )
            stub.setCompressionLevel(compressionLevel);
        if(indexOnTheFly != null)
            stub.setIndexOnTheFly(indexOnTheFly);
        if(generateMD5 != null)
            stub.setGenerateMD5(generateMD5);
        if(simplifyBAM != null)
            stub.setSimplifyBAM(simplifyBAM);

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

    private ArgumentDefinition disableWriteIndexArgumentDefinition(ArgumentSource source) {
        return new ArgumentDefinition( ArgumentIOType.ARGUMENT,
                                       boolean.class,
                                       DISABLE_INDEXING_FULLNAME,
                                       null,
                                       "Turn off on-the-fly creation of indices for output BAM files.",
                                       false,
                                       true,
                                       false,
                                       source.isHidden(),
                                       null,
                                       null,
                                       null,
                                       null );
    }

    private ArgumentDefinition enableMD5GenerationArgumentDefinition(ArgumentSource source) {
        return new ArgumentDefinition( ArgumentIOType.ARGUMENT,
                                       boolean.class,
                                       ENABLE_MD5_FULLNAME,
                                       null,
                                       "Enable on-the-fly creation of md5s for output BAM files.",
                                       false,
                                       true,
                                       false,
                                       source.isHidden(),
                                       null,
                                       null,
                                       null,
                                       null );
    }


    private ArgumentDefinition createSimplifyBAMArgumentDefinition(ArgumentSource source) {
        return new ArgumentDefinition( ArgumentIOType.ARGUMENT,
                                       boolean.class,
                                       SIMPLIFY_BAM_FULLNAME,
                                       SIMPLIFY_BAM_SHORTNAME,
                                       "If provided, output BAM files will be simplified to include just key reads for downstream variation discovery analyses (removing duplicates, PF-, non-primary reads), as well stripping all extended tags from the kept reads except the read group identifier",
                                       false,
                                       true,
                                       false,
                                       source.isHidden(),
                                       null,
                                       null,
                                       null,
                                       null );
    }
}
