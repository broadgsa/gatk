package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.calls.GenotypeCall;

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
     * Add a genotype, given a genotype locus
     *
     * @param locus the locus to add
     */
    @Override
    public void addGenotypeCall(GenotypeCall locus) {
        /*outStream.println(String.format("%s %s %c %s %s %f %f %f %f %d %s",
	                                        locus.getLocation(),
											"NOT OUTPUTED",
	                                        locus.getReferencebase(),
	                                        locus.getGenotypes().get(1).getBases().,
											genotype(),
	                                        qhat,
	                                        qstar,
                                            lodVsRef,
                                            lodVsNextBest,
	                                        depth,
											bases);   */
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
}
