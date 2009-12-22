package org.broadinstitute.sting.utils.genotype.glf;

import org.broadinstitute.sting.utils.genotype.GenotypeWriter;

/**
 * An extension of eth GenotypeWriter interface with support
 * for adding header lines.
 *
 * @author mhanna
 * @version 0.1
 */
public interface GLFGenotypeWriter extends GenotypeWriter {
    /**
     * Append the given header text to the GLF file.
     * @param headerText the file header to write out
     */    
    public void writeHeader(String headerText);

    /**
     * add a GLF record to the output file
     *
     * @param contigName   the contig name
     * @param contigLength the contig length
     * @param rec          the GLF record to write.
     */
    public void addGLFRecord(String contigName, int contigLength, GLFRecord rec);    
}
