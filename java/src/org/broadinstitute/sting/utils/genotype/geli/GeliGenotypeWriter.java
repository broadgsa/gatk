package org.broadinstitute.sting.utils.genotype.geli;

import org.broadinstitute.sting.utils.genotype.GenotypeWriter;
import net.sf.samtools.SAMFileHeader;
import edu.mit.broad.picard.genotype.geli.GenotypeLikelihoods;

/**
 * An extension of eth GenotypeWriter interface with support
 * for adding a header.
 *
 * @author mhanna
 * @version 0.1
 */
public interface GeliGenotypeWriter extends GenotypeWriter {
    /**
     * Write the file header.
     * @param fileHeader SAM file header from which to derive the geli header.
     */
    public void writeHeader(final SAMFileHeader fileHeader);

    /**
     * Writes the genotype likelihoods to the output.
     * @param gl genotype likelihoods to write.
     */
    public void addGenotypeLikelihoods(GenotypeLikelihoods gl);
}
