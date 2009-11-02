package org.broadinstitute.sting.utils.genotype.geli;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;


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
     * Add a genotype, given a call
     *
     * @param call the call to add
     */
    public void addGenotypeCall(Genotype call) {
        if ( !(call instanceof GeliGenotypeCall) )
            throw new IllegalArgumentException("Only GeliGenotypeCalls should be passed in to the Geli writers");
        GeliGenotypeCall gCall = (GeliGenotypeCall)call;

        char ref = gCall.getReference();

        double[] posteriors = gCall.getPosteriors();
        double[] lks;
        lks = Arrays.copyOf(posteriors, posteriors.length);
        Arrays.sort(lks);

        double nextVrsBest = lks[9] - lks[8];
        double nextVrsRef = 0;
        if (ref != 'X')
            nextVrsRef = lks[9] - posteriors[DiploidGenotype.createHomGenotype(ref).ordinal()];

        double maxMappingQual = 0;
        List<SAMRecord> recs = gCall.getReads();
        int readDepth = recs.size();
        for (SAMRecord rec : recs) {
            if (maxMappingQual < rec.getMappingQuality()) maxMappingQual = rec.getMappingQuality();
        }

        mWriter.println(String.format("%s    %16d  %c  %8d  %.0f  %s %.6f %.6f    %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f",
                                      gCall.getLocation().getContig(),
                                      gCall.getLocation().getStart(),
                                      ref,
                                      readDepth,
                                      maxMappingQual,
                                      gCall.getBases(),
                                      nextVrsRef,
                                      nextVrsBest,
                                      posteriors[0],
                                      posteriors[1],
                                      posteriors[2],
                                      posteriors[3],
                                      posteriors[4],
                                      posteriors[5],
                                      posteriors[6],
                                      posteriors[7],
                                      posteriors[8],
                                      posteriors[9]));
    }

    /**
     * add a no call to the genotype file, if supported.
     *
     * @param position the position to add the no call at
     */
    public void addNoCall(int position) {
        throw new UnsupportedOperationException("Geli text format doesn't support a no-call call.");
    }

    /** finish writing, closing any open files. */
    public void close() {
        mWriter.close();
    }

    /**
     * add a multi-sample call if we support it
     *
     * @param genotypes the list of genotypes
     */
    public void addMultiSampleCall(List<Genotype> genotypes, GenotypeLocusData metadata) {
        throw new UnsupportedOperationException("Geli text doesn't support multisample calls");
    }

    /** @return true if we support multisample, false otherwise */
    public boolean supportsMultiSample() {
        return false;
    }
}
