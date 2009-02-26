/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.sam;

import java.nio.ByteBuffer;

/**
 * Converter between binary and text CIGAR representation.
 */
class BinaryCigarCodec {
    private static final BinaryCigarCodec singleton = new BinaryCigarCodec();

    /**
     * It is not necssary to get the singleton but it is preferrable to use the same one
     * over and over vs. creating a new object for each BAMRecord.
     */
    static BinaryCigarCodec getSingleton() {
        return singleton;
    }

    int[] encode(final Cigar cigar) {
        if (cigar.numCigarElements() == 0) {
            return new int[0];
        }

        // Binary rep can be no longer than 1/2 of text rep
        // Although this is documented as uint, I think lengths will never get that long,
        // and it's a pain in Java.
        final int[] binaryCigar = new int[cigar.numCigarElements()];
        int binaryCigarLength = 0;
        for (int i = 0; i < cigar.numCigarElements(); ++i) {
            final CigarElement cigarElement = cigar.getCigarElement(i);
            final int op = CigarOperator.enumToBinary(cigarElement.getOperator());
            binaryCigar[binaryCigarLength++] = cigarElement.getLength() << 4 | op;
        }
        return binaryCigar;
    }

    Cigar decode(final ByteBuffer binaryCigar) {
        final Cigar ret = new Cigar();
        while (binaryCigar.hasRemaining()) {
            final int cigarette = binaryCigar.getInt();
            ret.add(binaryCigarToCigarElement(cigarette));
        }
        return ret;
    }

    Cigar decode(final int[] binaryCigar) {
        final Cigar ret = new Cigar();
        for (final int cigarette : binaryCigar) {
            ret.add(binaryCigarToCigarElement(cigarette));
        }
        return ret;
    }

    private static CigarElement binaryCigarToCigarElement(final int cigarette) {
        final int binaryOp = cigarette & 0xf;
        final int length = cigarette >> 4;
        return new CigarElement(length, CigarOperator.binaryToEnum(binaryOp));
    }
}
