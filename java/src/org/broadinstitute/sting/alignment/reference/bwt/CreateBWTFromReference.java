/*
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITHoc THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.alignment.reference.bwt;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.util.StringUtil;

import java.io.*;
import java.util.TreeSet;
import java.util.Comparator;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.alignment.reference.packing.PackUtils;

/**
 * Create a suffix array data structure.
 *
 * @author mhanna
 * @version 0.1
 */
public class CreateBWTFromReference {
    private String loadReference( File inputFile ) {
        // Read in the first sequence in the input file
        ReferenceSequenceFile reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(inputFile);
        ReferenceSequence sequence = reference.nextSequence();
        return StringUtil.bytesToString(sequence.getBases());
    }

    private String loadReverseReference( File inputFile ) {
        ReferenceSequenceFile reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(inputFile);
        ReferenceSequence sequence = reference.nextSequence();
        PackUtils.reverse(sequence.getBases());
        return StringUtil.bytesToString(sequence.getBases());
    }

    private Counts countOccurrences( String sequence ) {
        Counts occurrences = new Counts();
        for( char base: sequence.toCharArray() )
            occurrences.increment((byte)base);
        return occurrences;
    }

    private long[] createSuffixArray( String sequence ) {
        TreeSet<Integer> suffixArrayBuilder = new TreeSet<Integer>( new SuffixArrayComparator(sequence) );

        // Build out the suffix array using a custom comparator.
        System.out.printf("Creating sequence array of length %d%n", sequence.length() );
        for( int i = 0; i <= sequence.length(); i++ ) {
            suffixArrayBuilder.add(i);
            if( i % 100000 == 0 )
                System.out.printf("Added sequence %d%n", i);
        }

        // Copy the suffix array into an int array.
        long[] suffixArray = new long[suffixArrayBuilder.size()];
        int i = 0;
        for( Integer element: suffixArrayBuilder )
            suffixArray[i++] = element;

        return suffixArray;
    }

    private long[] invertSuffixArray( long[] suffixArray ) {
        long[] inverseSuffixArray = new long[suffixArray.length];
        for( int i = 0; i < suffixArray.length; i++ )
            inverseSuffixArray[(int)suffixArray[i]] = i;
        return inverseSuffixArray;
    }

    private long[] createCompressedSuffixArray( int[] suffixArray, int[] inverseSuffixArray ) {
        long[] compressedSuffixArray = new long[suffixArray.length];
        compressedSuffixArray[0] = inverseSuffixArray[0];
        for( int i = 1; i < suffixArray.length; i++ )
            compressedSuffixArray[i] = inverseSuffixArray[suffixArray[i]+1];
        return compressedSuffixArray;
    }

    private long[] createInversedCompressedSuffixArray( int[] compressedSuffixArray ) {
        long[] inverseCompressedSuffixArray = new long[compressedSuffixArray.length];
        for( int i = 0; i < compressedSuffixArray.length; i++ )
            inverseCompressedSuffixArray[compressedSuffixArray[i]] = i;
        return inverseCompressedSuffixArray;
    }

    private byte[] createBWT( String sequence, long[] suffixArray ) {
        byte[] bwt = new byte[suffixArray.length-1];
        int i = 0;
        for( long suffixArrayEntry: suffixArray ) {
            if( suffixArrayEntry == 0 )
                continue;
            bwt[i++] = (byte)sequence.charAt((int)suffixArrayEntry-1);
        }
        return bwt;
    }

    public static void main( String argv[] ) throws IOException {
        if( argv.length != 5 ) {
            System.out.println("USAGE: CreateBWTFromReference <input>.fasta <output bwt> <output rbwt> <output sa> <output rsa>");
            return;
        }

        String inputFileName = argv[0];
        File inputFile = new File(inputFileName);

        String bwtFileName = argv[1];
        File bwtFile = new File(bwtFileName);

        String rbwtFileName = argv[2];
        File rbwtFile = new File(rbwtFileName);

        String saFileName = argv[3];
        File saFile = new File(saFileName);

        String rsaFileName = argv[4];
        File rsaFile = new File(rsaFileName);

        CreateBWTFromReference creator = new CreateBWTFromReference();

        String sequence = creator.loadReference(inputFile);
        String reverseSequence = creator.loadReverseReference(inputFile);

        // Count the occurences of each given base.
        Counts occurrences = creator.countOccurrences(sequence);
        System.out.printf("Occurrences: a=%d, c=%d, g=%d, t=%d%n",occurrences.getCumulative(Bases.A),
                                                                  occurrences.getCumulative(Bases.C),
                                                                  occurrences.getCumulative(Bases.G),
                                                                  occurrences.getCumulative(Bases.T));

        // Generate the suffix array and print diagnostics.
        long[] suffixArrayData = creator.createSuffixArray(sequence);
        long[] reverseSuffixArrayData = creator.createSuffixArray(reverseSequence);

        // Invert the suffix array and print diagnostics.
        long[] inverseSuffixArray = creator.invertSuffixArray(suffixArrayData);
        long[] reverseInverseSuffixArray = creator.invertSuffixArray(reverseSuffixArrayData);

        SuffixArray suffixArray = new SuffixArray( inverseSuffixArray[0], occurrences, suffixArrayData );
        SuffixArray reverseSuffixArray = new SuffixArray( reverseInverseSuffixArray[0], occurrences, reverseSuffixArrayData );

        /*
        for( int i = 0; i < 8; i++ )
            System.out.printf("suffixArray[%d] = %d (%s...)%n", i, suffixArray.sequence[i], sequence.substring(suffixArray.sequence[i],Math.min(suffixArray.sequence[i]+100,sequence.length())));
        for( int i = 0; i < 8; i++ )
            System.out.printf("inverseSuffixArray[%d] = %d (%s...)%n", i, inverseSuffixArray[i], sequence.substring(i,Math.min(i+100,sequence.length())));
        */

        /*
        // Create the data structure for the compressed suffix array and print diagnostics.
        int[] compressedSuffixArray = creator.createCompressedSuffixArray(suffixArray.sequence,inverseSuffixArray);
        int reconstructedInverseSA = compressedSuffixArray[0];
        for( int i = 0; i < 8; i++ ) {
            System.out.printf("compressedSuffixArray[%d] = %d (SA-1[%d] = %d)%n", i, compressedSuffixArray[i], i, reconstructedInverseSA);
            reconstructedInverseSA = compressedSuffixArray[reconstructedInverseSA];
        }

        // Create the data structure for the inverse compressed suffix array and print diagnostics.
        int[] inverseCompressedSuffixArray = creator.createInversedCompressedSuffixArray(compressedSuffixArray);
        for( int i = 0; i < 8; i++ ) {
            System.out.printf("inverseCompressedSuffixArray[%d] = %d%n", i, inverseCompressedSuffixArray[i]);
        }
        */

        // Create the BWT.
        BWT bwt = new BWT(inverseSuffixArray[0], occurrences, creator.createBWT(sequence, suffixArray.sequence));
        BWT reverseBWT = new BWT( reverseInverseSuffixArray[0], occurrences, creator.createBWT(reverseSequence, reverseSuffixArray.sequence));

        byte[] bwtSequence = bwt.getSequence();
        System.out.printf("BWT: %s... (length = %d)%n", new String(bwtSequence,0,80),bwt.length());

        BWTWriter bwtWriter = new BWTWriter(bwtFile);
        bwtWriter.write(bwt);
        bwtWriter.close();

        BWTWriter reverseBWTWriter = new BWTWriter(rbwtFile);
        reverseBWTWriter.write(reverseBWT);
        reverseBWTWriter.close();

        /*
        SuffixArrayWriter saWriter = new SuffixArrayWriter(saFile);
        saWriter.write(suffixArray);
        saWriter.close();

        SuffixArrayWriter reverseSAWriter = new SuffixArrayWriter(rsaFile);
        reverseSAWriter.write(reverseSuffixArray);
        reverseSAWriter.close();
        */

        File existingBWTFile = new File(inputFileName+".bwt");
        BWTReader existingBWTReader = new BWTReader(existingBWTFile);
        BWT existingBWT = existingBWTReader.read();

        byte[] existingBWTSequence = existingBWT.getSequence();
        System.out.printf("Existing BWT: %s... (length = %d)%n",new String(existingBWTSequence,0,80),existingBWT.length());

        for( int i = 0; i < bwt.length(); i++ ) {
            if( bwtSequence[i] != existingBWTSequence[i] )
                throw new StingException("BWT mismatch at " + i);
        }

        File existingSAFile = new File(inputFileName+".sa");
        SuffixArrayReader existingSuffixArrayReader = new SuffixArrayReader(existingSAFile,existingBWT);
        SuffixArray existingSuffixArray = existingSuffixArrayReader.read();

        for(int i = 0; i < suffixArray.length(); i++) {
            if( i % 10000 == 0 )
                System.out.printf("Validating suffix array entry %d%n", i);
            if( suffixArray.get(i) != existingSuffixArray.get(i) )
                throw new StingException(String.format("Suffix array mismatch at %d; SA is %d; should be %d",i,existingSuffixArray.get(i),suffixArray.get(i)));
        }
    }

    /**
     * Compares two suffix arrays of the given sequence.  Will return whichever string appears
     * first in lexicographic order.
     */
    public static class SuffixArrayComparator implements Comparator<Integer> {
        /**
         * The data source for all suffix arrays.
         */
        private final String sequence;

        /**
         * Create a new comparator.
         * @param sequence Reference sequence to use as basis for comparison.
         */
        public SuffixArrayComparator( String sequence ) {
            this.sequence = sequence;
        }

        /**
         * Compare the two given suffix arrays.  Criteria for comparison is the lexicographic order of
         * the two substrings sequence[lhs:], sequence[rhs:].
         * @param lhs Left-hand side of comparison.
         * @param rhs Right-hand side of comparison.
         * @return How the suffix arrays represented by lhs, rhs compare.
         */
        public int compare( Integer lhs, Integer rhs ) {
            String lhsSuffixArray = sequence.substring(lhs);
            String rhsSuffixArray = sequence.substring(rhs);
            return lhsSuffixArray.compareTo(rhsSuffixArray);
        }
    }

}
