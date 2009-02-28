/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/

package edu.mit.broad.picard.importer.genotype;

import java.io.Closeable;
import java.io.File;

import edu.mit.broad.picard.PicardException;
import edu.mit.broad.sam.util.BinaryCodec;

/**
 *
 *
 * @author Doug Voet
 */
public class BedFileReader implements Closeable {
    private static final int LOWEST_2_BIT_MASK = 3; // binary 11
    private static final short BED_MAGIC_NUMBER = 7020;
//    private static final short BED_MAGIC_NUMBER = Short.parseShort("0110110000011011", 2);
    
    public static final byte MODE_INDIVIDUAL_MAJOR = 0;
    public static final byte MODE_SNP_MAJOR = 1;
    
    public static final byte GENOTYPE_AA = 0; // binary 00
    public static final byte GENOTYPE_NO_CALL = 1; // binary 01
    public static final byte GENOTYPE_AB = 2; // binary 10
    public static final byte GENOTYPE_BB = 3; // binary 11
    
    private final byte mode;
    private final BinaryCodec codec;
    private byte currentBlock;
    private int genotypeCount = 0;
    
    public BedFileReader(File bedFile) {
        this.codec = new BinaryCodec(bedFile, false);
        short fileMagicNumber = this.codec.readShort();
        if (fileMagicNumber != BED_MAGIC_NUMBER) {
            this.codec.close();
            throw new PicardException("Given file [" + bedFile.getAbsolutePath() + 
                    "] is not in bed file format... magic number does not match");
        }
        this.mode = codec.readByte();
    }
    
    public byte getMode() {
        return mode;    
    }

    @Override
    public void close() {
        this.codec.close();
    }
    
    public byte nextGenotype() {
        // there are 4 genotypes per byte so get a new byte every 4 genotypes read
        if (this.genotypeCount++ % 4 == 0) {
            this.currentBlock = this.codec.readByte();
        }
        
        // the 2 lowest order bits of currentBlock are the next genotype, pop them off
        byte genotype = (byte) (LOWEST_2_BIT_MASK & this.currentBlock);
        this.currentBlock >>>= 2;
        
        return genotype;
    }
    
    /**
     * Call this method when moving on to the next individual (in indiv-major mode) or next
     * snp (in snp-major mode).
     */
    public void dropRemainingBlock() {
        this.genotypeCount = 0;
    }
}
