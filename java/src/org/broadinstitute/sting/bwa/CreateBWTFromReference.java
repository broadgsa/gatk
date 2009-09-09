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

package org.broadinstitute.sting.bwa;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.util.StringUtil;

import java.io.*;
import java.util.TreeSet;
import java.util.Comparator;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * Create a suffix array data structure.
 *
 * @author mhanna
 * @version 0.1
 */
public class CreateBWTFromReference {
    private static final int ALPHABET_SIZE = 4;

    private String loadReference( File inputFile ) {
        // Read in the first sequence in the input file
        ReferenceSequenceFile reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(inputFile);
        ReferenceSequence sequence = reference.nextSequence();
        return StringUtil.bytesToString(sequence.getBases());
    }

    private int[] countOccurrences( String sequence ) {
        int occurrences[] = new int[ALPHABET_SIZE];
        for( char base: sequence.toCharArray() )
            occurrences[PackUtils.packBase((byte)base)]++;

        // Make occurrences cumulative
        for( int i = 1; i < ALPHABET_SIZE; i++ )
            occurrences[i] += occurrences[i-1];

        return occurrences;
    }

    private int[] createSuffixArray( String sequence ) {
        TreeSet<Integer> suffixArrayBuilder = new TreeSet<Integer>( new SuffixArrayComparator(sequence) );

        // Build out the suffix array using a custom comparator.
        System.out.printf("Creating sequence array of length %d%n", sequence.length() );
        for( int i = 0; i <= sequence.length(); i++ ) {
            suffixArrayBuilder.add(i);
            if( i % 100000 == 0 )
                System.out.printf("Added sequence %d%n", i);
        }

        // Copy the suffix array into an int array.
        int[] suffixArray = new int[suffixArrayBuilder.size()];
        int i = 0;
        for( Integer element: suffixArrayBuilder )
            suffixArray[i++] = element;

        return suffixArray;
    }

    private int[] invertSuffixArray( int[] suffixArray ) {
        int[] inverseSuffixArray = new int[suffixArray.length];
        for( int i = 0; i < suffixArray.length; i++ )
            inverseSuffixArray[suffixArray[i]] = i;
        return inverseSuffixArray;
    }

    private int[] createCompressedSuffixArray( int[] suffixArray, int[] inverseSuffixArray ) {
        int[] compressedSuffixArray = new int[suffixArray.length];
        compressedSuffixArray[0] = inverseSuffixArray[0];
        for( int i = 1; i < suffixArray.length; i++ )
            compressedSuffixArray[i] = inverseSuffixArray[suffixArray[i]+1];
        return compressedSuffixArray;
    }

    private int[] createInversedCompressedSuffixArray( int[] compressedSuffixArray ) {
        int[] inverseCompressedSuffixArray = new int[compressedSuffixArray.length];
        for( int i = 0; i < compressedSuffixArray.length; i++ )
            inverseCompressedSuffixArray[compressedSuffixArray[i]] = i;
        return inverseCompressedSuffixArray;
    }

    private byte[] createBWT( String sequence, int[] suffixArray ) {
        byte[] bwt = new byte[suffixArray.length-1];
        int i = 0;
        for( int suffixArrayEntry: suffixArray ) {
            if( suffixArrayEntry == 0 )
                continue;
            bwt[i++] = (byte)sequence.charAt(suffixArrayEntry-1);
        }
        return bwt;
    }

    public static void main( String argv[] ) throws IOException {
        if( argv.length != 3 ) {
            System.out.println("USAGE: CreateBWTFromReference <input>.fasta <output bwt> <output sa>");
            return;
        }

        String inputFileName = argv[0];
        File inputFile = new File(inputFileName);

        String bwtFileName = argv[1];
        File bwtFile = new File(bwtFileName);

        String saFileName = argv[2];
        File saFile = new File(saFileName);

        CreateBWTFromReference creator = new CreateBWTFromReference();

        String sequence = creator.loadReference(inputFile);

        // Generate the suffix array and print diagnostics.
        int[] suffixArray = creator.createSuffixArray(sequence);
        for( int i = 0; i < 8; i++ )
            System.out.printf("suffixArray[%d] = %d (%s...)%n", i, suffixArray[i], sequence.substring(suffixArray[i],Math.min(suffixArray[i]+100,sequence.length())));

        // Invert the suffix array and print diagnostics.
        int[] inverseSuffixArray = creator.invertSuffixArray(suffixArray);
        for( int i = 0; i < 8; i++ )
            System.out.printf("inverseSuffixArray[%d] = %d (%s...)%n", i, inverseSuffixArray[i], sequence.substring(i,Math.min(i+100,sequence.length())));

        // Create the data structure for the compressed suffix array and print diagnostics.
        int[] compressedSuffixArray = creator.createCompressedSuffixArray(suffixArray,inverseSuffixArray);
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

        // Count the occurences of each given base.
        int[] occurrences = creator.countOccurrences(sequence);
        System.out.printf("Occurrences: a=%d, c=%d, g=%d, t=%d%n",occurrences[0],occurrences[1],occurrences[2],occurrences[3]);

        // Create the BWT.
        byte[] bwt = creator.createBWT(sequence, suffixArray);

        String bwtAsString = new String(bwt);
        System.out.printf("BWT: %s...%n", bwtAsString.substring(0,80));

        OutputStream bwtOutputStream = new BufferedOutputStream(new FileOutputStream(bwtFile));

        ByteBuffer buffer = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
        buffer.putInt(inverseSuffixArray[0]);
        bwtOutputStream.write(buffer.array());
        bwtOutputStream.flush();

        PackedIntOutputStream occurrenceWriter = new PackedIntOutputStream(bwtOutputStream);
        occurrenceWriter.write(occurrences);
        occurrenceWriter.flush();

        BasePackedOutputStream<Integer> sequenceOutputStream = new BasePackedOutputStream<Integer>(Integer.class,bwtOutputStream,ByteOrder.LITTLE_ENDIAN);
        sequenceOutputStream.write(bwt);
        sequenceOutputStream.close();

        OutputStream saOutputStream = new BufferedOutputStream(new FileOutputStream(saFile));
        PackedIntOutputStream saIntWriter = new PackedIntOutputStream(saOutputStream);

        // SA file format is 'primary' (= SA-1[0]?), occurrence array, interval, sequence length, SA[]
        saIntWriter.write(inverseSuffixArray[0]);
        saIntWriter.write(occurrences);
        saIntWriter.write(1);
        saIntWriter.write(suffixArray.length-1);
        saIntWriter.write(suffixArray, 1, suffixArray.length-1);

        saIntWriter.close();

        File existingBwtFile = new File(inputFileName+".bwt");
        WordPackedInputStream inputStream = new WordPackedInputStream(existingBwtFile,ByteOrder.LITTLE_ENDIAN);
        byte[] existingBwt = inputStream.read();

        String existingBwtAsString = new String(existingBwt);
        System.out.printf("Existing BWT: %s...%n",existingBwtAsString.substring(0,80));

        for( int i = 0; i < bwt.length; i++ ) {
            if( bwt[i] != existingBwt[i] ) {
                System.out.printf("First bwt mismatch: %d%n",i);
                break;
            }
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
