package org.broadinstitute.sting.utils.genotype.geli;

import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidGenotype;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.genotype.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Arrays;


/**
 * @author aaron
 *         <p/>
 *         Class GeliTextWriter
 *         <p/>
 *         write out the geli text file format containing genotype information
 */
public class GeliTextWriter implements GenotypeWriter {
    // where we write to
    PrintWriter mWriter;

    /**
     * create a geli text writer
     *
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
    public void addGenotypeCall(Genotype locus) {
        double likelihoods[];
        int readDepth = -1;
        double nextVrsBest = 0;
        double nextVrsRef = 0;

        char ref = locus.getReference();

        if (locus instanceof ReadBacked) {
            readDepth = ((ReadBacked)locus).getReadCount();
        }
        if (!(locus instanceof GenotypesBacked)) {
            likelihoods = new double[10];
            Arrays.fill(likelihoods, 0.0);
        } else {
            likelihoods = ((LikelihoodsBacked) locus).getLikelihoods();
            double[] lks;
            lks = Arrays.copyOf(likelihoods,likelihoods.length);
            Arrays.sort(lks);
            nextVrsBest = lks[9] - lks[8];
            if (ref != 'X')  {
                int index = (DiploidGenotype.valueOf(Utils.dupString(ref,2)).ordinal());
                nextVrsRef = lks[9] - likelihoods[index];
            }
        }
        // we have to calcuate our own
        
        mWriter.println(String.format("%s    %16d  %c  %8d  %d  %s %.6f %.6f    %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f",
                                      locus.getLocation().getContig(),
                                      locus.getLocation().getStart(),
                                      ref,
                                      readDepth,
                                      -1,
                                      locus.getBases(),
                                      nextVrsRef,
                                      nextVrsBest,
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
