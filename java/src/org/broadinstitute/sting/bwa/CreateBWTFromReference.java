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

import java.io.IOException;
import java.io.File;
import java.util.TreeSet;
import java.util.Comparator;

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
            occurrences[ CreatePACFromReference.getPackedRepresentation((byte)base) ]++;

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

    private char[] createBWT( String sequence, int[] suffixArray ) {
        char[] bwt = new char[suffixArray.length];
        for( int i = 0; i < suffixArray.length; i++ ) {
            // Find the first character after the current character in the rotation.  If the character is past the end
            // (in other words, '$'), back up to the previous character.
            int sequenceEnd = Math.min((suffixArray[i]+suffixArray.length-1)%suffixArray.length, sequence.length()-1 );
            bwt[i] = sequence.charAt(sequenceEnd);
        }
        return bwt;
    }

    public static void main( String argv[] ) throws IOException {
        if( argv.length != 1 ) {
            System.out.println("No reference");
            return;
        }

        String inputFileName = argv[0];
        File inputFile = new File(inputFileName);

        CreateBWTFromReference creator = new CreateBWTFromReference();

        String sequence = creator.loadReference(inputFile);

        // Count the occurences of each given base.
        int[] occurrences = creator.countOccurrences(sequence);
        System.out.printf("Occurrences: a=%d, c=%d, g=%d, t=%d%n",occurrences[0],occurrences[1],occurrences[2],occurrences[3]);

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

        // Create the BWT.
        char[] bwt = creator.createBWT(sequence, suffixArray);
        System.out.printf("BWT: %s%n", new String(bwt));
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
