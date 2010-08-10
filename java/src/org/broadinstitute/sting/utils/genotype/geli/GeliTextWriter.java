package org.broadinstitute.sting.utils.genotype.geli;

import edu.mit.broad.picard.genotype.geli.GenotypeLikelihoods;
import net.sf.samtools.SAMFileHeader;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.StingException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.io.PrintWriter;


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

    // used to store the max mapping quality as a field in variant contexts
    public static final String MAXIMUM_MAPPING_QUALITY_ATTRIBUTE_KEY = "MAXIMUM_MAPPING_QUALITY";
    // used to store the max mapping quality as a field in variant contexts
    public static final String READ_COUNT_ATTRIBUTE_KEY = "READ_COUNT";

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
     * @param refBase required by the inteface; not used by this writer.
     */
    public void add(VariantContext vc, byte refBase) {
        throw new UnsupportedOperationException("We no longer support writing Geli");
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
