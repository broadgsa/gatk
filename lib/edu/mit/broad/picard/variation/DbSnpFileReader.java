/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2008 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.variation;

import java.io.*;
import java.util.*;
import edu.mit.broad.sam.SAMSequenceRecord;
import edu.mit.broad.sam.util.BinaryCodec;
import edu.mit.broad.picard.PicardException;
import edu.mit.broad.picard.io.IoUtil;

/**
 * Reader for DbSnp binary files.  See DbSnpFileGenerator for file format.
 */
public class DbSnpFileReader implements KnownVariantIterator
{
    private BinaryCodec codec = null;
    private KnownVariantCodec kvCodec = new KnownVariantCodec();
    List<SAMSequenceRecord> dictionary; 
    private Map<Integer,SAMSequenceRecord> refIndexToName = new HashMap<Integer,SAMSequenceRecord>();
    private KnownVariant next = null;
    private int dbSnpCount = -1;

    /**
     * Constructor
     *
     * @param dbSnpFile  The binary dbSnp file to read
     */
    public DbSnpFileReader(File dbSnpFile)
    {
        codec = new BinaryCodec(new DataInputStream(IoUtil.openFileForReading(dbSnpFile)));
        readHeader();
        next = readNextDbSnp();
    }

    /**
     * Returns an iterator over a set of elements of type KnownVariant.
     *
     * @return an Iterator
     */
    public Iterator<KnownVariant> iterator()
    {
        return this;
    }

    /**
     * Returns true if the iteration has more elements.
     *
     * @return  true if the iterator has more elements.
     */
    public boolean hasNext()
    {
        return next != null;
    }

    /**
     * Returns the next element in the iteration.
     *
     * @return the next KnownVariant in the iteratoion
     */
    public KnownVariant next()
    {
        if (!hasNext()) throw new NoSuchElementException();
        KnownVariant result = next;
        next = readNextDbSnp();
        return result;
    }

    /** Allows peeking at the next value without advaning the iterator. */
    public KnownVariant peek() {
        return this.next;
    }

    /**
     * Not supported.
     *
     * @throws UnsupportedOperationException    
     */
    public void remove()
    {
        throw new UnsupportedOperationException("Remove() not supported.");
    }

    /**
     * Closes the underlying stream, via the BinaryCodec's close() method
     */
    public void close()
    {
            codec.close();
    }

    /**
     * Reads the header data from the binary file, validates the version, and populates <code>refIndexToName</code>
     * 
     * @throws IOException
     */
    private void readHeader()
    {
        // Verify that we are using the correct version
        String ver = kvCodec.decodeMagicNumber(codec);
        if (!ver.equals(KnownVariantCodec.MAGIC_NUMBER))
        {
            throw new RuntimeException("Unsupported dbSnp file version: " + ver);
        }

        // Read the number of reference sequences and then the sequences themselves
        dictionary = kvCodec.decodeSequenceDictionary(codec);
        for (int i = 0; i < dictionary.size(); i++)
        {
            refIndexToName.put(i, dictionary.get(i));
        }

        dbSnpCount = codec.readInt();
    }

    /**
     * Reads the next dbSnp record from the binary file
     *
     * @return  the populated KnownVariant object
     */
    private KnownVariant readNextDbSnp() {
        KnownVariant kv = kvCodec.decodeKnownVariant(codec);
        if (kv != null) {
            kv.setRefrenceSequence(refIndexToName.get(kv.getSequenceIndex()).getSequenceName());
        }
        return kv;
    }

    /**
     * Returns the SequenceDictionary for this file in SAM format
     *
     * @return an ordered List of SAMSequenceRecords
     */
    public List<SAMSequenceRecord> getSequenceDictionary() { return dictionary; }

    /**
     * Returns the total number of dbSnp records encoded in the file
     *
     * @return  total dbSnps encoded in the file
     */
    public int getCountDbSnpRecords() { return dbSnpCount; }
}