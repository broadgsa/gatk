package org.broadinstitute.sting.utils.genotype;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.geli.*;
import org.broadinstitute.sting.utils.genotype.glf.*;
import org.broadinstitute.sting.utils.genotype.vcf.*;

import java.io.File;
import java.io.PrintStream;
import java.util.Set;


/**
 * @author aaron
 *         <p/>
 *         Class GenotypeWriterFactory
 *         <p/>
 *         A descriptions should go here. Blame aaron if it's missing.
 */
public class GenotypeWriterFactory {
    /** available genotype writers */
    public enum GENOTYPE_FORMAT {
        GELI, GLF, GFF, TABULAR, GELI_BINARY, VCF
    }

    /**
     * create a genotype writer
     * @param format the format
     * @param destination the destination file
     * @return the genotype writer object
     */
    public static GenotypeWriter create(GENOTYPE_FORMAT format, File destination) {
        switch (format) {
            case GLF:
                return new GLFWriter(destination);
            case GELI:
                return new GeliTextWriter(destination);
            case GELI_BINARY:
                return new GeliAdapter(destination);
            case VCF:
                return new VCFGenotypeWriterAdapter(destination);
            default:
                throw new StingException("Genotype writer " + format.toString() + " is not implemented");
        }
    }

    public static GenotypeWriter create(GENOTYPE_FORMAT format, PrintStream destination) {
        switch (format) {
            case GELI:
                return new GeliTextWriter(destination);
            case GLF:
                return new GLFWriter(destination);
            case VCF:
                return new VCFGenotypeWriterAdapter(destination);
            default:
                throw new StingException("Genotype writer to " + format.toString() + " to standard output is not implemented");
        }
    }

    public static void writeHeader(GenotypeWriter writer,
                                   SAMFileHeader header,
                                   Set<String> sampleNames,
                                   Set<VCFHeaderLine> headerInfo) {
        // VCF
        if ( writer instanceof VCFGenotypeWriter ) {
            ((VCFGenotypeWriter)writer).writeHeader(sampleNames, headerInfo);
        }
        // GELI 
        else if ( writer instanceof GeliGenotypeWriter ) {
            ((GeliGenotypeWriter)writer).writeHeader(header);
        }
        // GLF
        else if ( writer instanceof GLFGenotypeWriter ) {
            ((GLFGenotypeWriter)writer).writeHeader(header.toString());
        }        
        // nothing to do for GELI TEXT
    }

    /**
     * create a genotype call
     * @param format the format
     * @param ref the reference base
     * @param loc the location
     * @return an unpopulated genotype call object
     */
    public static GenotypeCall createSupportedGenotypeCall(GENOTYPE_FORMAT format, char ref, GenomeLoc loc) {
        switch (format) {
            case VCF:
                return new VCFGenotypeCall(ref, loc);
            case GELI:
            case GELI_BINARY:
                return new GeliGenotypeCall(ref, loc);
            case GLF:
                return new GLFGenotypeCall(ref, loc);
            default:
                throw new StingException("Genotype format " + format.toString() + " is not implemented");
        }
    }

    /**
     * create a genotype locus data object
     * @param format the format
     * @param ref the reference base
     * @param loc the location
     * @param type the variant type
     * @return an unpopulated genotype locus data object
     */
    public static VariationCall createSupportedCall(GENOTYPE_FORMAT format, char ref, GenomeLoc loc, Variation.VARIANT_TYPE type) {
        switch (format) {
            case VCF:
                return new VCFVariationCall(ref, loc, type);
            case GELI:
            case GELI_BINARY:
                return null;
            case GLF:
                return null;
            default:
                throw new StingException("Genotype format " + format.toString() + " is not implemented");
        }
    }
}
