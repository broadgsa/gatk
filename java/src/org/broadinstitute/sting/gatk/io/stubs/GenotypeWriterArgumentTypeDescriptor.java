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
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;
import org.broadinstitute.sting.utils.genotype.GenotypeWriterFactory;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;

import java.io.File;
import java.util.List;
import java.util.Arrays;

/**
 * Injects new command-line arguments into the system providing support for the genotype writer.
 *
 * @author mhanna
 * @version 0.1
 */
public class GenotypeWriterArgumentTypeDescriptor extends ArgumentTypeDescriptor {
    /**
     * The engine into which output stubs should be fed.
     */
    private GenomeAnalysisEngine engine;

    /**
     * Create a new GenotypeWriter argument, notifying the given engine when that argument has been created.
     * @param engine the engine to be notified.
     */
    public GenotypeWriterArgumentTypeDescriptor(GenomeAnalysisEngine engine) {
        this.engine = engine;
    }

    /**
     * Reports whether this ArgumentTypeDescriptor supports the given type.
     * @param type The type to check.
     * @return True if the argument is a GenotypeWriter.
     */
    @Override
    public boolean supports( Class type ) {
        return GenotypeWriter.class.equals(type);
    }

    /**
     * Create the argument definitions associated with this source.
     * Assumes that this type descriptor is relevant for this source.
     * @param source Source class and field for the given argument.
     * @return A list of all associated argument definitions.
     */
    @Override
    public List<ArgumentDefinition> createArgumentDefinitions( ArgumentSource source ) {
        return Arrays.asList( createGenotypeFileArgumentDefinition(source),
                              createGenotypeFormatArgumentDefinition(source) );
    }

    /**
     * This command-line argument descriptor does want to override the provided default value.
     * @return true always.
     */
    @Override
    public boolean overridesDefault() {
        return true;
    }

    /**
     * Provide the default value for this argument.
     * @return A VCFGenotypeWriter which writes to the default output stream.
     */
    @Override
    public Object getDefault() {
        GenotypeWriterStub defaultGenotypeWriter = new VCFGenotypeWriterStub(engine,System.out);
        engine.addOutput(defaultGenotypeWriter);
        return defaultGenotypeWriter;       
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
        String writerFileName = getArgumentValue(createGenotypeFileArgumentDefinition(source),matches);
        File writerFile = writerFileName != null ? new File(writerFileName) : null;

        // Get the format of the genotype object, if it exists.
        String genotypeFormatText = getArgumentValue(createGenotypeFormatArgumentDefinition(source),matches);
        GenotypeWriterFactory.GENOTYPE_FORMAT genotypeFormat = GenotypeWriterFactory.GENOTYPE_FORMAT.VCF;
        if(genotypeFormatText != null) {
            try {
                genotypeFormat = Enum.valueOf(GenotypeWriterFactory.GENOTYPE_FORMAT.class,genotypeFormatText);
            }
            catch(IllegalArgumentException ex) {
                throw new StingException(String.format("Genotype format %s is invalid.",genotypeFormatText));
            }            
        }

        // Create a stub for the given object.
        GenotypeWriterStub stub;
        switch(genotypeFormat) {
            case GELI:
                stub = (writerFile != null) ? new GeliTextGenotypeWriterStub(engine, writerFile) : new GeliTextGenotypeWriterStub(engine,System.out);
                break;
            case GELI_BINARY:
                if(writerFile == null)
                    throw new StingException("Geli binary files cannot be output to the console.");
                stub = new GeliBinaryGenotypeWriterStub(engine, writerFile);
                break;
            case GLF:
                stub = (writerFile != null) ? new GLFGenotypeWriterStub(engine, writerFile) : new GLFGenotypeWriterStub(engine,System.out);
                break;
            case VCF:
                stub = (writerFile != null) ? new VCFGenotypeWriterStub(engine, writerFile) : new VCFGenotypeWriterStub(engine,System.out);
                break;
            default:
                throw new StingException("Unable to create stub for file format " + genotypeFormat);
        }

        engine.addOutput(stub);

        return stub;
    }

    /**
     * Gets the definition of the argument representing the BAM file itself.
     * @param source Argument source for the BAM file.  Must not be null.
     * @return Argument definition for the BAM file itself.  Will not be null.
     */
    private ArgumentDefinition createGenotypeFileArgumentDefinition(ArgumentSource source) {
        ArgumentDescription description = this.getArgumentDescription(source);

        boolean isFullNameProvided = description.fullName().trim().length() > 0;
        boolean isShortNameProvided = description.shortName().trim().length() > 0;

        String fullName = isFullNameProvided ? description.fullName().trim() : "variants_out";

        // If the short name is provided, use that.  If the user hasn't provided any names at all, use
        // the default.  If somewhere in the middle, leave the short name blank.
        String shortName;
        if( isShortNameProvided )
            shortName = description.shortName().trim();
        else if( !isFullNameProvided )
            shortName = "varout";
        else
            shortName = null;

        return new ArgumentDefinition( getIOType(source),
                                       fullName,
                                       shortName,
                                       getDoc(source),
                                       isRequired(source),
                                       false,
                                       source.isMultiValued(),
                                       getExclusiveOf(source),
                                       getValidationRegex(source),
                                       null );
    }

    /**
     * Creates the optional compression level argument for the BAM file.
     * @param source Argument source for the BAM file.  Must not be null.
     * @return Argument definition for the BAM file itself.  Will not be null.
     */
    private ArgumentDefinition createGenotypeFormatArgumentDefinition(ArgumentSource source) {
        return new ArgumentDefinition( getIOType(source),
                                       "variant_output_format",
                                       "vf",
                                       "Format to be used to represent variants; default is VCF",
                                       false,
                                       false,
                                       false,
                                       null,
                                       null,
                                       null );
    }

}
