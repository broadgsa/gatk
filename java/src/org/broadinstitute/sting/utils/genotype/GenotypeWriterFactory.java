package org.broadinstitute.sting.utils.genotype;

import org.broad.tribble.vcf.VCFHeader;
import org.broadinstitute.sting.utils.vcf.GATKVCFWriter;
import org.broadinstitute.sting.utils.genotype.vcf.*;

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
        GELI, GLF, GELI_BINARY, VCF
    }

    /**
     * create a genotype writer
     * @param destination the destination file
     * @return the genotype writer object
     */
    public static GenotypeWriter create(File destination) {
        return new GATKVCFWriter(destination);
    }

    public static GenotypeWriter create(PrintStream destination) {
        return new GATKVCFWriter(destination);
    }

    public static void writeHeader(GenotypeWriter writer, VCFHeader vcfHeader) {
        ((VCFGenotypeWriter)writer).writeHeader(vcfHeader);
    }
}
