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

import java.io.*;
import java.util.*;
import edu.mit.broad.sam.SAMSequenceRecord;
import edu.mit.broad.sam.SAMFileReader;
import edu.mit.broad.sam.util.BinaryCodec;
import edu.mit.broad.picard.io.IoUtil;
import edu.mit.broad.picard.util.TabbedTextFileParser;
import edu.mit.broad.picard.util.Log;

/**
 * Generates a binary version of the data for all dbSnps from a UCSU snp###.txt file.  Files with SNP data
 * can be downloaded here:  http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/.  See KnownVariantCodec.java
 * for binary file format.
 */
public class DbSnpFileGenerator {
    // Codes from the DbSnp file that we will handle.  All others are ignored.
    // Package visibility for testing purposes.
    static final String snp       = "single";    // code in DbSnp file for a SNP
    static final String insertion = "insertion"; // code in DbSnp file for an insertion
    static final String deletion  = "deletion";  // code in DbSnp file for a deletion
    static final String indel     = "in-del";    // code in DbSnp file for an insertion/deletion

    private File snpFile;
    private File seqDictionaryFile;
    private Map<String, Integer> sequenceToIndex = new HashMap<String, Integer>();
    private List<SAMSequenceRecord> dictionary;
    private BinaryCodec codec;
    private KnownVariantCodec kvCodec = new KnownVariantCodec();
    private Map<String, SortedSet<KnownVariant>> sequenceToSnps;

    private final Log log = Log.getInstance(DbSnpFileGenerator.class);

    /**
     * Protected constructor so we can use a temporary file during testing
     * @param snpFile               The UCSC dbSnp file
     * @param seqDictionaryFile     The Sequence Dictionary
     * @param tempOutputFile            The binary file to write to
     */
    DbSnpFileGenerator(File snpFile, File seqDictionaryFile, File tempOutputFile) {
        this.snpFile = snpFile;
        this.seqDictionaryFile = seqDictionaryFile;
        this.codec = new BinaryCodec(new DataOutputStream(IoUtil.openFileForWriting(tempOutputFile)));
    }

    /**
     * Writes the full binary dbSnp file and calls close on the BinaryCodec.
     */
    public void writeDbSnpFile() {
        kvCodec.encode(KnownVariantCodec.MAGIC_NUMBER, codec);
        writeReferenceSequences();
        writeDbSnpRecords();
        codec.close();
    }

    /**
     * Writes the number of reference sequences and then the sequences themselves
     */
    private void writeReferenceSequences()  {
        SAMFileReader sam = new SAMFileReader(this.seqDictionaryFile);
        this.dictionary = sam.getFileHeader().getSequences();
        kvCodec.encode(this.dictionary, codec);
    }

    /**
     * Writes all the dbSnp records to the file in the order of the reference sequences
     * in the sequence dictionary file.
     */
    private void writeDbSnpRecords() {
        sequenceToSnps = new HashMap<String, SortedSet<KnownVariant>>();
        int count = 0;

        TabbedTextFileParser parser = new TabbedTextFileParser(true, snpFile);
        while(parser.hasNext())  {
            String parts[] = parser.next();
            String sequence = parts[1];

            // If we don't have this sequence in our dictionary, ignore it
            if (!getSequenceToIndex().containsKey(sequence)) {
                continue;
            }

            int start = Integer.parseInt(parts[2]) + 1; // We go from a zero-based to a 1-based system.
            int end = Integer.parseInt(parts[3]);

            String var = parts[11];

            // We only care about SNPs, insertions, and deletions; otherwise skip it
            VariantType type = null;
            if (var.equals(snp)) {
                type = VariantType.SNP;
                end = start;    // For SNPs, we mark the start and end as the same location
            }
            // For insertions and deletions, we mark the base on either side of the affected reference sequence
            else if (var.equals(insertion)) {
                type = VariantType.insertion;
                end = start + 1;    // Insertions are always length 1
            }
            else if (var.equals(deletion)) {
                type = VariantType.deletion;
                start = start - 1;
                end++;
            }
            else if (var.equals(indel)) {   // For indels, we do one each of an insertion (here) and a deletion (below)
                type = VariantType.insertion;
                start = start - 1;
                end = start + 1;
            }
            else {
                continue;
            }

            if (!sequenceToSnps.containsKey(sequence)) {
                sequenceToSnps.put(sequence, new TreeSet<KnownVariant>());
            }
            SortedSet<KnownVariant> sequenceVars = sequenceToSnps.get(sequence);

            boolean validated = !parts[12].equals("unknown");
            String name = parts[4];

            sequenceVars.add(new KnownVariant(name, getSequenceToIndex().get(sequence), start, end, type, validated));
            count++;

            // If it's an in-del, we add it as a deletion (in addition to the insertion we also added) so we
            // will have two records in our binary format for the one record in the text file
            if (var.equals(indel)) {
                sequenceVars.add(new KnownVariant(name, getSequenceToIndex().get(sequence), start,
                        Integer.parseInt(parts[3])+1, VariantType.deletion, validated));
                count++;
            }
        }

        codec.writeInt(count);
        // Loop through the sequences from the sequence dictionary in order
        for (int i = 0; i < dictionary.size(); i++) {
            // And write their known variants in order
            if (sequenceToSnps.containsKey(dictionary.get(i).getSequenceName())) {
                for (Iterator<KnownVariant> it = sequenceToSnps.get(dictionary.get(i).getSequenceName()).iterator();
                     it.hasNext();) {
                    kvCodec.encode(it.next(), codec);
                }
            }
        }
        log.info("Wrote " + count + " dbSnp records.");
    }

    /**
     * Returns the map of sequences to their index in the reference dictionary,
     * creating it if it does not already exist
     *
     * @return the map of sequences to their index in the reference dictionary
     */
    private Map<String, Integer> getSequenceToIndex() {
        if (sequenceToIndex.keySet().size() == 0) {
            for (int i = 0; i < dictionary.size(); i++) {
                sequenceToIndex.put(dictionary.get(i).getSequenceName(), i);
            }
        }
        return sequenceToIndex;
    }

}