package org.broadinstitute.sting.utils.genotype;

import net.sf.samtools.SAMFileHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFGenotypeRecord;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.geli.*;
import org.broadinstitute.sting.utils.genotype.glf.*;
import org.broadinstitute.sting.utils.genotype.vcf.*;

import java.io.File;
import java.io.PrintStream;
import java.util.Set;
import java.util.HashSet;


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
        GELI, GLF, GELI_BINARY, VCF
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
            if ( headerInfo == null )
                headerInfo = new HashSet<VCFHeaderLine>(VCFGenotypeRecord.getSupportedHeaderStrings());
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
}
