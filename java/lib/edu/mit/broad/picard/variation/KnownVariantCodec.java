/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.variation;

import edu.mit.broad.sam.SAMSequenceRecord;
import edu.mit.broad.sam.util.BinaryCodec;
import edu.mit.broad.sam.util.RuntimeEOFException;

import java.util.ArrayList;
import java.util.List;

/**
 * Class for encoding and deconding binary data about KnownVariants
 *
 * IMPORTANT!  This class assumes that a KnownVariant instance is 1-based and end-inclusive
 * and that the binary format is 0-based and end-exclusive.
 *
 * The format for the binary dbSnp file is as follows:
 *
 *  Field        Description                                             Type            Value
 *  -----        -----------                                             ----            -----
 *  magic        Known variant magic number                              char[4]         DBS\1
 *  n_ref        # reference sequences                                   int32
 *
 *      -- List of references information (n = n_ref)
 *    l_name     length of the reference name plus 1 (including NULL)    int32
 *    name       Name; NULL terminated                                   char[l_name]
 *    t_ref      Length of the reference sequence                        int32
 *
 *
 *  n_snps       # of Known Variant records                              int32
 *
 *      -- List of DBSnps
 *    block_size Length of the remainder of the block
 *    rID        Reference sequence ID (-1 <= rId <= n_ref)              int32
 *    pos        0-based leftmost coordinate                             int32
 *    snp_len    Length of the dbSnp                                     int32
 *    type       type of SNP                                             int8            0 = deletion
 *                                                                                       1 = het
 *                                                                                       2 = in-del
 *                                                                                       3 = insertion
 *                                                                                       4 = microsatellite
 *                                                                                       5 = mixed
 *                                                                                       6 = mnp
 *                                                                                       7 = named
 *                                                                                       8 = single
 *                                                                                       9 = unknown
 *    validated  whether the SNP has been validated                      int8            1 | 0
 *    name       name of the dbSnp; NULL terminated                      char[block_size-15]
 *
 *  @author Kathleen Tibbetts
 **/
public class KnownVariantCodec
{
    public static final String MAGIC_NUMBER = "DBS\1";
    private static final int KV_RECORD_LENGTH_LESS_NAME = 15;
    
    /**
     * Reads data about a known variant from the BinaryCodec and instantiates a KnownVariant
     * object with those values
     *
     * @param codec     The BinaryCodec from which to read
     * @return a populated KnownVariant object
     */
    public KnownVariant decodeKnownVariant(BinaryCodec codec)
    {
        int blockSize;
        try {
            blockSize = codec.readInt();
        }
        catch (RuntimeEOFException e) {
            return null;
        }
        int seqIndex = codec.readInt();
        int startPos = codec.readInt() + 1; // Switch to 1-based
        int endPos = codec.readInt();
        byte[] buffer = new byte[1];
        codec.readBytes(buffer);
        VariantType type = VariantType.getVariantTypeFromOrdinal((int) buffer[0]);
        codec.readBytes(buffer);
        boolean validated = ((int) buffer[0]) == 1;
        String name = codec.readString(blockSize - KV_RECORD_LENGTH_LESS_NAME);
        codec.readBytes(buffer); // Skip the null terminator
        return new KnownVariant(name, seqIndex, startPos, endPos, type, validated);

    }

    /**
     * Writes data from a KnownVariant in the expected format to the BinaryCodec
     * 
     * @param variant  The KnownVariant to encode
     * @param codec    The BinaryCodec to which to write
     */
    public void encode(KnownVariant variant, BinaryCodec codec)
    {
        codec.writeInt(variant.getName().length() + KV_RECORD_LENGTH_LESS_NAME);// Length of the rest of the block
        codec.writeInt(variant.getSequenceIndex());              // Index of the reference sequence
        codec.writeInt((int)variant.getStartPos()-1);            // Switch to 0-based leftmost coordinate
        codec.writeInt((int)variant.getEndPos());                // end position, exclusive
        byte b[] = new byte[1];
        b[0] = (byte)variant.getType().ordinal();                // Type
        codec.writeBytes(b);
        b[0] = (byte)(variant.isValidated() ? 1 : 0);            // Validated
        codec.writeBytes(b);
        codec.writeString(variant.getName(), false, true);       // The null-terminated name
    }

    /**
     * Reads data about the Sequence Dictionary from the BinaryCodec and instantiates a List of
     * SAMSequenceRecords with those values
     *
     * @param codec     The BinaryCodec from which to read
     * @return a populated List of SAMSequenceRecords
     */
    public List<SAMSequenceRecord> decodeSequenceDictionary(BinaryCodec codec)
    {
        int total = codec.readInt();
        List<SAMSequenceRecord> dictionary = new ArrayList<SAMSequenceRecord>(total);
        for (int i = 0; i < total; i++)
        {
            int len = codec.readInt();
            // Read the name, leaving off and then skipping the null terminator
            String name = codec.readString(len-1);
            byte[] buffer = new byte[1];
            codec.readBytes(buffer);
            int seqLength = codec.readInt();
            SAMSequenceRecord rec = new SAMSequenceRecord(name);
            rec.setSequenceLength(seqLength);
            dictionary.add(rec);
        }
        return dictionary;
    }

    /**
     * Writes a Sequence Dictionary in the format excpected to the BinaryCodec
     *
     * @param dictionary  The list of SAMSequenceRecords to encode
     * @param codec    The BinaryCodec to which to write
     */
    public void encode(List<SAMSequenceRecord> dictionary, BinaryCodec codec)
    {
        codec.writeInt(dictionary.size());
        for (SAMSequenceRecord sequence : dictionary)
        {
            codec.writeString(sequence.getSequenceName(), true, true);
            codec.writeInt(sequence.getSequenceLength());
        }

    }

    /**
     * Reads data about the Magic Number from the BinaryCodec and returns a string with its value
     *
     * @param codec     The BinaryCodec from which to read
     * @return a Magic Number
     */
    public String decodeMagicNumber(BinaryCodec codec)
    {
        return codec.readString(4);
    }

    /**
     * Writes a Magic Number in the format excpected to the BinaryCodec
     *
     * @param magicNumber  The magic number to encode
     * @param codec    The BinaryCodec to which to write
     */
    public void encode(String magicNumber, BinaryCodec codec)
    {
        codec.writeString(magicNumber, false, false);
    }
}
