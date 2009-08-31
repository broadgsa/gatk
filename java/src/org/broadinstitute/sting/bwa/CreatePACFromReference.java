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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.bwa;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.picard.reference.ReferenceSequence;

import java.io.*;

import org.broadinstitute.sting.utils.StingException;

/**
 * Generate a .PAC file from a given reference.
 *
 * @author hanna
 * @version 0.1
 */

public class CreatePACFromReference {
    private static final int ALPHABET_SIZE = 4;

    public static void main( String argv[] ) throws FileNotFoundException, IOException {
        if( argv.length != 1 ) {
            System.out.println("No reference");
            return;
        }

        // Read in the first sequence in the input file
        String inputFileName = argv[0];
        File inputFile = new File(inputFileName);
        ReferenceSequenceFile reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(inputFile);
        ReferenceSequence sequence = reference.nextSequence();

        // Target file for output
        File outputFile = new File(inputFileName + ".mypac");
        BufferedOutputStream outputStream = new BufferedOutputStream(new FileOutputStream(outputFile));

        // Number of bytes, rounded up, plus one extra byte indicating how many bases the last byte holds.
        byte packed = 0;
        int positionInPack = 0;

        for( byte base: sequence.getBases() ) {
            // Pack base into the appropriate bits of the byte.
            packed |= (getPackedRepresentation(base) << 2*(ALPHABET_SIZE-positionInPack-1));

            // Increment the packed counter.  If all possible bases have been squeezed into this byte, write it out.
            positionInPack = ++positionInPack % 4;
            if( positionInPack == 0 ) {
                outputStream.write(packed);
                packed = 0;
            }
        }

        // Last (incomplete) block in file.
        if( positionInPack > 0 )
            outputStream.write(packed);

        // Last character of a .pac file is how many bases are in the last byte.
        outputStream.write(positionInPack);
        
        outputStream.close();
    }

    /**
     * Gets the two-bit representation of a base.  A=00b, C=01b, G=10b, T=11b.
     * @param base ASCII value for the base to pack.
     * @return A byte from 0-3 indicating the base's packed value.
     */
    private static byte getPackedRepresentation(byte base) {
        switch( base ) {
            case 'A':
                return 0;
            case 'C':
                return 1;
            case 'G':
                return 2;
            case 'T':
                return 3;
            default:
                throw new StingException("Unknown base type: " + base);
        }
    }
}
