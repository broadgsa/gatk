package org.broadinstitute.sting.utils.genotype;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.utils.genotype.glf.GLFWriter;

import java.io.File;


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
        GELI, GLF, GFF, TABULAR;
    }

    /**
     * create a genotype writer
     * @param format the format
     * @param header the sam file header
     * @param destination the destination file
     * @return the genotype writer object
     */
    public static GenotypeWriter create(GENOTYPE_FORMAT format, SAMFileHeader header, File destination) {
        switch (format) {
            case GLF:
                return new GLFWriter(header.toString(), destination);
            case GELI:
                return new GeliAdapter(destination, header);
        }
        return null;
    }
}
