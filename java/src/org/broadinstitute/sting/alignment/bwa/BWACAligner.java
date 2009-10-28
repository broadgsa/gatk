package org.broadinstitute.sting.alignment.bwa;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.alignment.Alignment;

import java.io.File;

/**
 * An aligner using the BWA/C implementation.
 *
 * @author mhanna
 * @version 0.1
 */
public class BWACAligner {
    static {
        System.loadLibrary("bwa");
    }

    /**
     * A pointer to the C++ object representing the BWA engine.
     */
    private long thunkPointer = 0;

    /**
     * Create a pointer to the BWA/C thunk.
     * @param annFileName Name of the ann file?
     * @param ambFileName Name of the amb file?
     * @param pacFileName Packed representation of the forward reference.
     * @param forwardBWTFileName Name of the file where the forward BWT is stored.
     * @param forwardSAFileName Name of te file where the forward suffix array is stored.
     * @param reverseBWTFileName Name of the file where the reverse BWT is stored.
     * @param reverseSAFileName Name of the file where the reverse SA is stored.
     * @return Pointer to the BWA/C thunk.
     */
    protected native long create(String annFileName,
                                 String ambFileName,
                                 String pacFileName,
                                 String forwardBWTFileName,
                                 String forwardSAFileName,
                                 String reverseBWTFileName,
                                 String reverseSAFileName);

    /**
     * Destroy the 
     * @param thunkPointer Pointer to the allocated thunk.
     */
    protected native void destroy(long thunkPointer);

    public BWACAligner(String annFileName,
                     String ambFileName,
                     String pacFileName,
                     String forwardBWTFileName,
                     String forwardSAFileName,
                     String reverseBWTFileName,
                     String reverseSAFileName) {
        if(thunkPointer != 0)
            throw new StingException("BWA/C attempting to reinitialize.");
        thunkPointer = create(annFileName,
                              ambFileName,
                              pacFileName,
                              forwardBWTFileName,
                              forwardSAFileName,
                              reverseBWTFileName,
                              reverseSAFileName);
    }

    /**
     * Close this instance of the BWA pointer and delete its resources.
     */
    public void close() {
        if(thunkPointer == 0)
            throw new StingException("BWA/C close attempted, but BWA/C was never properly initialized.");
        destroy(thunkPointer);
    }

    /**
     * Align the given base array to the BWT.  The base array should be in ASCII format.
     * @param bases ASCII representation of byte array.
     * @return an array of indices into the bwa.
     */
    public Alignment[] align(byte[] bases) {
        if(thunkPointer == 0)
            throw new StingException("BWA/C align attempted, but BWA/C was never properly initialized.");
        return align(thunkPointer,bases);
    }

    /**
     * Caller to the .  The base array should be in
     * ASCII format.
     * @param thunkPointer pointer to the C++ object managing BWA/C.
     * @param bases ASCII representation of byte array.
     */
    protected native Alignment[] align(long thunkPointer, byte[] bases);

    public static void main(String[] args) {
        String prefix = "/Users/mhanna/reference/Ecoli/Escherichia_coli_K12_MG1655.fasta";
        BWACAligner thunk = new BWACAligner(prefix + ".ann",
                                        prefix + ".amb",
                                        prefix + ".pac",
                                        prefix + ".bwt",
                                        prefix + ".sa",
                                        prefix + ".rbwt",
                                        prefix + ".rsa");

        SAMFileReader reader = new SAMFileReader(new File("/Users/mhanna/reference/Ecoli/MV1994.bam"));
        reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        int count = 0;

        for(SAMRecord read: reader) {
            count++;
            //if(count > 1) break;
            Alignment[] alignments = thunk.align(read.getReadBases());
            System.out.printf("Read: %s: ", read.getReadName());
            for(Alignment alignment: alignments)
                System.out.printf("tig = %d; pos = %d, neg strand = %b, mapQ = %d, cigar = %s;",
                                  alignment.getContigIndex(),
                                  alignment.getAlignmentStart(),
                                  alignment.isNegativeStrand(),
                                  alignment.getMappingQuality(),
                                  alignment.getCigarString());
            if(count % 10000 == 0) System.out.printf("Processed %d reads.%n",count);
            System.out.printf("%n");
        }

        thunk.close();
    }
}
