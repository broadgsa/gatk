package org.broadinstitute.sting.utils.genotype.gff;

import org.broadinstitute.sting.utils.genotype.GenotypeWriter;
import org.broadinstitute.sting.utils.genotype.LikelihoodObject;
import org.broadinstitute.sting.utils.genotype.IndelLikelihood;
import net.sf.samtools.SAMSequenceRecord;


/**
 * @author aaron
 *         <p/>
 *         Class GFFWriter
 *         <p/>
 *         A descriptions should go here. Blame aaron if it's missing.
 */
public class GFFWriter implements GenotypeWriter {

    /**
     * add a single point genotype call to the file
     *
     * @param contig        the contig you're calling in
     * @param position      the position on the contig
     * @param referenceBase the reference base
     * @param readDepth     the read depth at the specified position
     * @param likelihoods   the likelihoods of each of the possible alleles
     */
    @Override
    public void addGenotypeCall(SAMSequenceRecord contig, int position, float rmsMapQuals, char referenceBase, int readDepth, LikelihoodObject likelihoods) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * add a variable length call to the genotype file
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
        //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * add a no call to the genotype file, if supported.
     *
     * @param position  the position
     * @param readDepth the read depth
     */
    @Override
    public void addNoCall(int position, int readDepth) {
        //To change body of implemented methods use File | Settings | File Templates.
    }


    /** finish writing, closing any open files. */
    @Override
    public void close() {
        //To change body of implemented methods use File | Settings | File Templates.
    }
}
