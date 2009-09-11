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
import java.nio.ByteOrder;

/**
 * Generate a .PAC file from a given reference.
 *
 * @author hanna
 * @version 0.1
 */

public class CreatePACFromReference {
    public static void main( String argv[] ) throws IOException {
        if( argv.length != 3 ) {
            System.out.println("USAGE: CreatePACFromReference <input>.fasta <output pac> <output rpac>");
            return;
        }

        // Read in the first sequence in the input file
        String inputFileName = argv[0];
        File inputFile = new File(inputFileName);
        ReferenceSequenceFile reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(inputFile);
        ReferenceSequence sequence = reference.nextSequence();

        // Target file for output
        writeSequence( new File(argv[1]), sequence.getBases() );

        // Reverse the bases in the reference
        PackUtils.reverse(sequence.getBases());

        // Target file for output
        writeSequence( new File(argv[2]), sequence.getBases() );
    }

    private static void writeSequence( File outputFile, byte[] bases ) throws IOException {
        OutputStream outputStream = new FileOutputStream(outputFile);

        BasePackedOutputStream<Byte> basePackedOutputStream = new BasePackedOutputStream<Byte>(Byte.class, outputStream, ByteOrder.BIG_ENDIAN);
        basePackedOutputStream.write(bases);

        outputStream.write(bases.length%PackUtils.ALPHABET_SIZE);

        outputStream.close();
    }
}
