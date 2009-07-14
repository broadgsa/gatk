package org.broadinstitute.sting.utils.genotype;

import net.sf.samtools.SAMSequenceRecord;
import org.broadinstitute.sting.utils.StingException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;


/**
 * 
 * @author aaron 
 * 
 * Class TabularLF
 *
 * the tabular likelihood format, as an implementation of the genotype interface
 */
public class TabularLFWriter implements GenotypeWriter {
    /**
     * where to print the tabular genotype likelihood info to
     */
    public PrintStream outStream;

    /**
     * construct, writing to a specified file
     * @param writeTo
     */
    public TabularLFWriter(File writeTo) {
        try {
            outStream = new PrintStream(writeTo);
        } catch (FileNotFoundException e) {
            throw new StingException("Unable to write to specified file: " + writeTo.getName());
        }
        // print the header out
        outStream.println("location sample_name ref alt genotype qhat qstar lodVsRef lodVsNextBest depth bases");
    }

    /**
     * add a single point genotype call to the
     *
     * @param contig        the contig you're calling in
     * @param position      the position on the contig
     * @param referenceBase the reference base
     * @param readDepth     the read depth at the specified position
     * @param likelihoods   the likelihoods of each of the possible alleles
     */
    @Override
    public void addGenotypeCall(SAMSequenceRecord contig, int position, float rmsMapQuals, char referenceBase, int readDepth, LikelihoodObject likelihoods) {
         /**return String.format("%s %s %c %c %s %f %f %f %f %d %s",
	                                        location,
											contig.getSpecies(),
	                                        ref,
	                                        alt,
											genotype(),
	                                        qhat,
	                                        qstar,
                                            lodVsRef,
                                            lodVsNextBest,
	                                        depth,
											bases);*/        
    }

    /**
     * add a variable length call to the genotyper
     *
     * @param contig        the contig you're calling in
     * @param position      the position on the genome
     * @param rmsMapQuals   the root mean square of the mapping qualities
     * @param readDepth     the read depth
     * @param refBase       the reference base
     * @param firstHomZyg   the first homozygous indel
     * @param secondHomZyg  the second homozygous indel (if present, null if not)
     * @param hetLikelihood the heterozygous likelihood
     */
    @Override
    public void addVariableLengthCall(SAMSequenceRecord contig, int position, float rmsMapQuals, int readDepth, char refBase, IndelLikelihood firstHomZyg, IndelLikelihood secondHomZyg, byte hetLikelihood) {
        throw new StingException("TabularLFWriter doesn't support variable length calls");
    }

    /**
     * add a no call to the genotype file, if supported.
     *
     * @param position
     * @param readDepth
     */
    @Override
    public void addNoCall(int position, int readDepth) {
        throw new StingException("TabularLFWriter doesn't support no-calls");
    }

    /** finish writing, closing any open files. */
    @Override
    public void close() {
        if (this.outStream != null) {
            outStream.close();
        }
    }
}
