package org.broadinstitute.sting.utils.genotype;

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
     * Add a genotype, given a genotype locus
     *
     * @param locus the locus to add
     */
    @Override
    public void addGenotypeCall(GenotypeOutput locus) {
        /**
         * This output is not correct, but I don't we even use this format anymore.  If we do, someone
         * should change this code
          */
        outStream.println(String.format("%s %s %c %s %s %f %f %f %f %d %s",
	                                        locus.getLocation().toString(),
											"NOT OUTPUTED",
	                                        locus.getReferencebase(),
	                                        locus.getBases(),
											locus.getBases(),
	                                        -1,
	                                        -1,
                                            locus.getBestRef(),
                                            locus.getBestNext(),
	                                        locus.getReadDepth(),
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
}
