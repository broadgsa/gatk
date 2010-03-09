package org.broadinstitute.sting.utils.genotype.tabular;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.genotype.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.Arrays;


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
     * @param writeTo file to write to
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
     * Add a genotype, given a variant context
     *
     * @param vc  the variant context representing the call to add
     */
    public void addCall(VariantContext vc) {
        if ( vc.getNSamples() != 1 )
            throw new IllegalArgumentException("The tabular LF format does not support multi-sample or no-calls");

        org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype genotype = vc.getGenotypes().values().iterator().next();
        if ( genotype.isNoCall() )
            throw new IllegalArgumentException("The tabular LF format does not support no-calls");

        ReadBackedPileup pileup;
        double[] likelihoods;
        if ( genotype instanceof CalledGenotype) {
            pileup = ((CalledGenotype)genotype).getReadBackedPileup();
            likelihoods = ((CalledGenotype)genotype).getLikelihoods();
        } else {
            pileup = (ReadBackedPileup)genotype.getAttribute(CalledGenotype.READBACKEDPILEUP_ATTRIBUTE_KEY);
            likelihoods = (double[])genotype.getAttribute(CalledGenotype.LIKELIHOODS_ATTRIBUTE_KEY);
        }

        if ( likelihoods == null ) {
            likelihoods = new double[10];
            Arrays.fill(likelihoods, Double.MIN_VALUE);
        }

        int readDepth = pileup == null ? -1 : pileup.size();

        double nextVrsBest = 0;
        double nextVrsRef = 0;
        char ref = vc.getReference().toString().charAt(0);


        /**
         * This output is not correct, but I don't we even use this format anymore.  If we do, someone
         * should change this code
          */
        outStream.println(String.format("%s %s %c %s %s %f %f %f %f %d %s",
	                                        vc.getLocation().toString(),
											"NOT OUTPUTED",
	                                        ref,
	                                        genotype.getGenotypeString(),
											genotype.getGenotypeString(),
	                                        -1,
	                                        -1,
                                            nextVrsRef,
                                            nextVrsBest,
	                                        readDepth,
											genotype.getGenotypeString()));
    }

    /** finish writing, closing any open files. */
    public void close() {
        if (this.outStream != null) {
            outStream.close();
        }
    }
}
