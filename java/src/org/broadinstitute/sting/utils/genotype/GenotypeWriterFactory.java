package org.broadinstitute.sting.utils.genotype;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.geli.GeliAdapter;
import org.broadinstitute.sting.utils.genotype.geli.GeliTextWriter;
import org.broadinstitute.sting.utils.genotype.glf.GLFWriter;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeWriterAdapter;

import java.io.File;
import java.io.PrintStream;


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
        GELI, GLF, GFF, TABULAR, GELI_BINARY, VCF;
    }

    /**
     * create a genotype writer
     * @param format the format
     * @param header the sam file header
     * @param destination the destination file
     * @return the genotype writer object
     */
    public static GenotypeWriter create(GENOTYPE_FORMAT format, SAMFileHeader header, File destination, String source, String referenceName ) {
        switch (format) {
            case GLF:
                return new GLFWriter(header.toString(), destination);
            case GELI:
                return new GeliTextWriter(destination);
            case GELI_BINARY:
                return new GeliAdapter(destination, header);
            case VCF:
                return new VCFGenotypeWriterAdapter(source, referenceName, destination);
            default:
                throw new StingException("Genotype writer " + format.toString() + " is not implemented");
        }
    }

    public static GenotypeWriter create(GENOTYPE_FORMAT format, SAMFileHeader header, PrintStream destination, String source, String referenceName ) {
        switch (format) {
            case GELI:
                return new GeliTextWriter(destination);
            case GLF:
                return new GLFWriter(header.toString(), destination);
            case VCF:
                return new VCFGenotypeWriterAdapter(source, referenceName, destination);
            default:
                throw new StingException("Genotype writer to " + format.toString() + " to standard output is not implemented");
        }
    }
}
