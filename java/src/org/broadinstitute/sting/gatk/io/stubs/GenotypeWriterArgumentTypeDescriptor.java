package org.broadinstitute.sting.gatk.io.stubs;

import org.broadinstitute.sting.utils.cmdLine.ArgumentTypeDescriptor;
import org.broadinstitute.sting.utils.cmdLine.ArgumentSource;
import org.broadinstitute.sting.utils.cmdLine.ArgumentMatches;
import org.broadinstitute.sting.utils.cmdLine.ArgumentDefinition;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;
import org.broadinstitute.sting.utils.genotype.GenotypeWriterFactory;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;

import java.io.File;
import java.util.List;
import java.util.Arrays;

import net.sf.samtools.SAMFileReader;

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
     * @param engine
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
     * Convert the given argument matches into a single object suitable for feeding into the ArgumentSource.
     * @param source Source for this argument.
     * @param type
     * @param matches Matches that match with this argument.
     * @return Transform from the matches into the associated argument.
     */
    @Override
    public Object parse( ArgumentSource source, Class type, ArgumentMatches matches )  {
        String writerFileName = getArgumentValue(createGenotypeFileArgumentDefinition(source),matches);
        if(writerFileName == null)
            throw new StingException("Genotype format was supplied, but no file was supplied to contain the genotype info..");

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

        GenotypeWriterStub stub = new GenotypeWriterStub(engine, new File(writerFileName),genotypeFormat);

        engine.addOutput(stub);

        return stub;
    }

    /**
     * Gets the definition of the argument representing the BAM file itself.
     * @param source Argument source for the BAM file.  Must not be null.
     * @return Argument definition for the BAM file itself.  Will not be null.
     */
    private ArgumentDefinition createGenotypeFileArgumentDefinition(ArgumentSource source) {
        Argument description = this.getArgumentDescription(source);

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
    private ArgumentDefinition createGenotypeFormatArgumentDefinition(ArgumentSource source) {
        return new ArgumentDefinition( "variant_output_format",
                                       "vf",
                                       "Format to be used to represent variants; default is VCF",
                                       false,
                                       false,
                                       false,
                                       null,
                                       null );
    }

}
