package org.broadinstitute.sting.utils.genotype.geli;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.GenotypeCall;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.List;


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
     * @param file
     */
    public GeliTextWriter(File file) {
        try {
            mWriter = new PrintWriter(file);
        } catch (FileNotFoundException e) {
            throw new StingException("Unable to open file " + file.toURI());
        }
        mWriter.println("#Sequence       Position        ReferenceBase   NumberOfReads   MaxMappingQuality       BestGenotype    BtrLod  BtnbLod dbSNP   AA      AC      AG      AT      CC      CG      CT      GG      GT      TT");
    }

    /**
     * Add a genotype, given a genotype locus
     *
     * @param locus the locus to add
     */
    public void addGenotypeCall(GenotypeCall locus) {
        if (locus.getPosteriors().size() != 10) throw new IllegalArgumentException("Geli text only supports SNP calls, with a diploid organism (i.e. posterior array size of 10)");


        // this is to perserve the format string that we used to use
        double[] likelihoods = new double[10];
        int index = 0;
        List<Genotype> lt = locus.getLexigraphicallySortedGenotypes();
        for (Genotype G: lt) {
            likelihoods[index] = G.getLikelihood();
            index++;
        }

        mWriter.println( String.format("%s    %16d  %c  %8d  %d  %s %.6f %.6f    %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f",
	                                        locus.getLocation().getContig(),
                                            locus.getLocation().getStart(),
											locus.getReferencebase(),
                                            locus.getReadDepth(),
                                            -1,
	                                        locus.getGenotypes().get(0).getBases(),
	                                        locus.getBestVrsRef().second.getScore(),
	                                        locus.getBestVrsNext().second.getScore(),
                                            likelihoods[0],
                                            likelihoods[1],
                                            likelihoods[2],
                                            likelihoods[3],
                                            likelihoods[4],
                                            likelihoods[5],
                                            likelihoods[6],
                                            likelihoods[7],
                                            likelihoods[8],
                                            likelihoods[9]));
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
