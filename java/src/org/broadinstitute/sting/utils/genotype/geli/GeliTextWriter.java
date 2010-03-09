package org.broadinstitute.sting.utils.genotype.geli;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.genotype.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collections;

import edu.mit.broad.picard.genotype.geli.GenotypeLikelihoods;


/**
 * @author aaron
 *         <p/>
 *         Class GeliTextWriter
 *         <p/>
 *         write out the geli text file format containing genotype information
 */
public class GeliTextWriter implements GeliGenotypeWriter {
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
    }

    public GeliTextWriter(PrintStream out) {
        mWriter = new PrintWriter(out);
    }

    public final static String headerLine = "#Sequence       Position        ReferenceBase   NumberOfReads   MaxMappingQuality       BestGenotype    BtrLod  BtnbLod    AA      AC      AG      AT      CC      CG      CT      GG      GT      TT";

    /**
     * Write the file header.
     * @param fileHeader SAM file header from which to derive the geli header.
     */
    public void writeHeader(final SAMFileHeader fileHeader) {
        // ignore the SAM header; the geli text header is fixed.
        mWriter.println(headerLine);        
        mWriter.flush();  // necessary so that writing to an output stream will work
    }

    /**
     * Add a genotype, given a variant context
     *
     * @param vc  the variant context representing the call to add
     */
    public void addCall(VariantContext vc) {

        char ref = vc.getReference().toString().charAt(0);

        if ( vc.getNSamples() != 1 )
            throw new IllegalArgumentException("The Geli format does not support multi-sample or no-calls");

        org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype genotype = vc.getGenotypes().values().iterator().next();
        if ( genotype.isNoCall() )
            throw new IllegalArgumentException("The Geli format does not support no-calls");

        ReadBackedPileup pileup;
        double[] posteriors;
        if ( genotype instanceof CalledGenotype ) {
            pileup = ((CalledGenotype)genotype).getReadBackedPileup();
            posteriors = ((CalledGenotype)genotype).getPosteriors();
        } else {
            pileup = (ReadBackedPileup)genotype.getAttribute(CalledGenotype.READBACKEDPILEUP_ATTRIBUTE_KEY);
            posteriors = (double[])genotype.getAttribute(CalledGenotype.POSTERIORS_ATTRIBUTE_KEY);
        }

        if ( posteriors == null )
            throw new IllegalArgumentException("The Geli format requires posteriors");

        double[] lks;
        lks = Arrays.copyOf(posteriors, posteriors.length);
        Arrays.sort(lks);

        double nextVrsBest = lks[9] - lks[8];
        double nextVrsRef = 0;
        if (ref != 'X')
            nextVrsRef = lks[9] - posteriors[DiploidGenotype.createHomGenotype(ref).ordinal()];

        int readCount = 0;
        double maxMappingQual = 0;
        if ( pileup != null ) {
            readCount = pileup.size();
            for (PileupElement p : pileup ) {
                if ( maxMappingQual < p.getMappingQual() )
                    maxMappingQual = p.getMappingQual();
            }
        }

        ArrayList<Character> alleles = new ArrayList<Character>();
        for ( Allele a : genotype.getAlleles() )
            alleles.add(a.toString().charAt(0));
        Collections.sort(alleles);
        StringBuffer sb = new StringBuffer();
        for ( Character base : alleles )
            sb.append(base);

        mWriter.println(String.format("%s    %16d  %c  %8d  %.0f  %s %.6f %.6f    %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f",
                                      vc.getLocation().getContig(),
                                      vc.getLocation().getStart(),
                                      ref,
                                      readCount,
                                      maxMappingQual,
                                      sb.toString(),
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
        mWriter.flush();  // necessary so that writing to an output stream will work
    }

    public void addGenotypeLikelihoods(GenotypeLikelihoods gl) {
        mWriter.println(gl.toString());
        mWriter.flush();  // necessary so that writing to an output stream will work
    }

    /** finish writing, closing any open files. */
    public void close() {
        mWriter.flush();
        mWriter.close();
    }
}
