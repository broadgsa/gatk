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

package org.broadinstitute.gatk.utils.clipping;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.CigarUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.Assert;

import java.util.LinkedList;
import java.util.List;
import java.util.Stack;

public class ReadClipperTestUtils {
    //Should contain all the utils needed for tests to mass produce
    //reads, cigars, and other needed classes

    final static byte [] BASES = {'A', 'C', 'T', 'G'};
    final static byte [] QUALS = {2, 15, 25, 30};
    final static String CIGAR = "4M";
    final static CigarElement[] cigarElements = { new CigarElement(1, CigarOperator.HARD_CLIP),
                                                  new CigarElement(1, CigarOperator.SOFT_CLIP),
                                                  new CigarElement(1, CigarOperator.INSERTION),
                                                  new CigarElement(1, CigarOperator.DELETION),
                                                  new CigarElement(1, CigarOperator.MATCH_OR_MISMATCH)};

    /**
     * Make a read fom the CIGAR
     *
     * @param cigar the CIGAR
     * @param lengthChange change in read length relative the CIGAR length
     * @return artificial read
     */
    public static GATKSAMRecord makeReadFromCigar(Cigar cigar, int lengthChange) {
        int readLength = cigar.getReadLength();
        if ( readLength >= -lengthChange ) {
            readLength += lengthChange;
        }

        return ArtificialSAMUtils.createArtificialRead(Utils.arrayFromArrayWithLength(BASES, readLength), Utils.arrayFromArrayWithLength(QUALS, readLength), cigar.toString());
    }

    /**
     * Make a read from the CIGAR string
     *
     * @param cigarString string used to create a CIGAR
     * @param lengthChange change in read length relative the CIGAR length
     * @return artificial read
     */
    public static GATKSAMRecord makeReadFromCigar(String cigarString, int lengthChange) {
        return makeReadFromCigar(CigarUtils.cigarFromString(cigarString), lengthChange);
    }

    public static List<Cigar> generateCigarList(int maximumLength) {
        return generateCigarList(maximumLength, cigarElements);
    }

        /**
        * This function generates every valid permutation of cigar strings (with a given set of cigarElement) with a given length.
        *
        * A valid cigar object obeys the following rules:
        *  - No Hard/Soft clips in the middle of the read
        *  - No deletions in the beginning / end of the read
        *  - No repeated adjacent element (e.g. 1M2M -> this should be 3M)
        *  - No consecutive I/D elements
        *
        * @param maximumLength the maximum number of elements in the cigar
        * @return a list with all valid Cigar objects
        */
    public static List<Cigar> generateCigarList(int maximumLength, CigarElement[] cigarElements) {
        int numCigarElements = cigarElements.length;
        LinkedList<Cigar> cigarList = new LinkedList<Cigar>();
        byte [] cigarCombination = new byte[maximumLength];

        Utils.fillArrayWithByte(cigarCombination, (byte) 0);               // we start off with all 0's in the combination array.
        int currentIndex = 0;
        while (true) {
            Cigar cigar = createCigarFromCombination(cigarCombination, cigarElements);    // create the cigar
            cigar = CigarUtils.combineAdjacentCigarElements(cigar);                   // combine adjacent elements
            if (CigarUtils.isCigarValid(cigar)) {                                     // check if it's valid
                cigarList.add(cigar);                                      // add it
            }

            boolean currentIndexChanged = false;
            while (currentIndex < maximumLength && cigarCombination[currentIndex] == numCigarElements - 1) {
                currentIndex++;                                            // find the next index to increment
                currentIndexChanged = true;                                // keep track of the fact that we have changed indices!
            }

            if (currentIndex == maximumLength)                             // if we hit the end of the array, we're done.
                break;

            cigarCombination[currentIndex]++;                              // otherwise advance the current index

            if (currentIndexChanged) {                                     // if we have changed index, then...
                for (int i = 0; i < currentIndex; i++)
                    cigarCombination[i] = 0;                               // reset everything from 0->currentIndex
                currentIndex = 0;                                          // go back to the first index
            }
        }

        return cigarList;
    }

    private static Cigar createCigarFromCombination(byte[] cigarCombination, CigarElement[] cigarElements) {
        Cigar cigar = new Cigar();
        for (byte i : cigarCombination) {
            cigar.add(cigarElements[i]);
        }
        return cigar;
    }

    public static GATKSAMRecord makeRead() {
        return ArtificialSAMUtils.createArtificialRead(BASES, QUALS, CIGAR);
    }

    /**
     * Asserts that the two reads have the same bases, qualities and cigar strings
     *
     * @param actual the calculated read
     * @param expected the expected read
     */
    public static void assertEqualReads(GATKSAMRecord actual, GATKSAMRecord expected) {
        // If they're both not empty, test their contents
        if(!actual.isEmpty() && !expected.isEmpty()) {
            Assert.assertEquals(actual.getReadBases(), expected.getReadBases());
            Assert.assertEquals(actual.getBaseQualities(), expected.getBaseQualities());
            Assert.assertEquals(actual.getCigarString(), expected.getCigarString());
        }
        // Otherwise test if they're both empty
        else
            Assert.assertEquals(actual.isEmpty(), expected.isEmpty());
     }
}
