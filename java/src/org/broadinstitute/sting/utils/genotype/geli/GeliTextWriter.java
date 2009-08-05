package org.broadinstitute.sting.utils.genotype.geli;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.GenotypeOutput;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;
import org.broadinstitute.sting.gatk.walkers.genotyper.SSGGenotypeCall;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.PrintStream;


/**
 * 
 * @author aaron 
 * 
 * Class GeliTextWriter
 *
 * A descriptions should go here. Blame aaron if it's missing.
 */
public class GeliTextWriter implements GenotypeWriter {
    // where we write to
    PrintWriter mWriter;

    /**
     * create a geli text writer
     * @param file the file to write to
     */
    public GeliTextWriter(File file) {
        try {
            mWriter = new PrintWriter(file);
        } catch (FileNotFoundException e) {
            throw new StingException("Unable to open file " + file.toURI());
        }
        mWriter.println(headerLine);
    }

    public GeliTextWriter(PrintStream out) {
        mWriter = new PrintWriter(out);
        mWriter.println(headerLine);
    }

    public final static String headerLine = "#Sequence       Position        ReferenceBase   NumberOfReads   MaxMappingQuality       BestGenotype    BtrLod  BtnbLod    AA      AC      AG      AT      CC      CG      CT      GG      GT      TT";

    /**
     * Add a genotype, given a genotype locus
     *
     * @param locus the locus to add
     */
    public void addGenotypeCall(GenotypeOutput locus) {
        SSGGenotypeCall call = (SSGGenotypeCall)locus;

        mWriter.println( String.format("%s    %16d  %c  %8d  %d  %s %.6f %.6f    %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f",
	                                        locus.getLocation().getContig(),
                                            locus.getLocation().getStart(),
											locus.getReferencebase(),
                                            call.getReadDepth(),
                                            -1,
	                                        locus.getBases(),
	                                        call.getBestRef(),
	                                        call.getBestNext(),
                                            call.getLikelihoods()[0],
                                            call.getLikelihoods()[1],
                                            call.getLikelihoods()[2],
                                            call.getLikelihoods()[3],
                                            call.getLikelihoods()[4],
                                            call.getLikelihoods()[5],
                                            call.getLikelihoods()[6],
                                            call.getLikelihoods()[7],
                                            call.getLikelihoods()[8],
                                            call.getLikelihoods()[9]));
    }

    /**
     * add a no call to the genotype file, if supported.
     *
     * @param position the position to add the no call at
     */
    @Override
    public void addNoCall(int position) {
        throw new UnsupportedOperationException("Geli text format doesn't support a no-call call.");
    }

    /** finish writing, closing any open files. */
    @Override
    public void close() {
        mWriter.close();
    }
}
