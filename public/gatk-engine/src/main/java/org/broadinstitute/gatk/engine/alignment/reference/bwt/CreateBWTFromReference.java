/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.alignment.reference.bwt;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.broadinstitute.gatk.engine.alignment.reference.packing.PackUtils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.io.File;
import java.io.IOException;

/**
 * Create a suffix array data structure.
 *
 * @author mhanna
 * @version 0.1
 */
public class CreateBWTFromReference {
    private byte[] loadReference( File inputFile ) {
        // Read in the first sequence in the input file
        ReferenceSequenceFile reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(inputFile);
        ReferenceSequence sequence = reference.nextSequence();
        return sequence.getBases();
    }

    private byte[] loadReverseReference( File inputFile ) {
        ReferenceSequenceFile reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(inputFile);
        ReferenceSequence sequence = reference.nextSequence();
        PackUtils.reverse(sequence.getBases());
        return sequence.getBases();
    }

    private Counts countOccurrences( byte[] sequence ) {
        Counts occurrences = new Counts();
        for( byte base: sequence )
            occurrences.increment(base);
        return occurrences;
    }

    private long[] createSuffixArray( byte[] sequence ) {
        return SuffixArray.createFromReferenceSequence(sequence).sequence;
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

        byte[] sequence = creator.loadReference(inputFile);
        byte[] reverseSequence = creator.loadReverseReference(inputFile);

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
        BWT bwt = BWT.createFromReferenceSequence(sequence);
        BWT reverseBWT = BWT.createFromReferenceSequence(reverseSequence);

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
                throw new ReviewedGATKException("BWT mismatch at " + i);
        }

        File existingSAFile = new File(inputFileName+".sa");
        SuffixArrayReader existingSuffixArrayReader = new SuffixArrayReader(existingSAFile,existingBWT);
        SuffixArray existingSuffixArray = existingSuffixArrayReader.read();

        for(int i = 0; i < suffixArray.length(); i++) {
            if( i % 10000 == 0 )
                System.out.printf("Validating suffix array entry %d%n", i);
            if( suffixArray.get(i) != existingSuffixArray.get(i) )
                throw new ReviewedGATKException(String.format("Suffix array mismatch at %d; SA is %d; should be %d",i,existingSuffixArray.get(i),suffixArray.get(i)));
        }
    }

}
