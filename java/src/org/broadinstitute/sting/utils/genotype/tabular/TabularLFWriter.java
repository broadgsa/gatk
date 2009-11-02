package org.broadinstitute.sting.utils.genotype.tabular;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;


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
     * Add a genotype, given a genotype locus
     *
     * @param locus the locus to add
     */
    @Override
    public void addGenotypeCall(Genotype locus) {
        double likelihoods[];
        int readDepth = -1;
        double nextVrsBest = 0;
        double nextVrsRef = 0;
        if (!(locus instanceof LikelihoodsBacked)) {
            likelihoods = new double[10];
            Arrays.fill(likelihoods, Double.MIN_VALUE);
        } else {
            likelihoods = ((LikelihoodsBacked) locus).getLikelihoods();

        }
        char ref = locus.getReference();

        if (locus instanceof ReadBacked) {
            readDepth = ((ReadBacked)locus).getReadCount();
        }
        /**
         * This output is not correct, but I don't we even use this format anymore.  If we do, someone
         * should change this code
          */
        outStream.println(String.format("%s %s %c %s %s %f %f %f %f %d %s",
	                                        locus.getLocation().toString(),
											"NOT OUTPUTED",
	                                        ref,
	                                        locus.getBases(),
											locus.getBases(),
	                                        -1,
	                                        -1,
                                            nextVrsRef,
                                            nextVrsBest,
	                                        readDepth,
											locus.getBases()));
    }

    /**
     * add a no call to the genotype file, if supported.
     *
     * @param position
     */
    @Override
    public void addNoCall(int position) {
        throw new StingException("TabularLFWriter doesn't support no-calls");
    }

    /** finish writing, closing any open files. */
    @Override
    public void close() {
        if (this.outStream != null) {
            outStream.close();
        }
    }


    /**
     * add a multi-sample call if we support it
     *
     * @param genotypes the list of genotypes, that are backed by sample information
     */
    @Override
    public void addMultiSampleCall(List<Genotype> genotypes, GenotypeLocusData metadata) {
        throw new UnsupportedOperationException("Tabular LF doesn't support multisample calls");
    }

    /** @return true if we support multisample, false otherwise */
    @Override
    public boolean supportsMultiSample() {
        return false;
    }
}
