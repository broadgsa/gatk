/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils.clipping;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.TextCigarCodec;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
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


    public static GATKSAMRecord makeReadFromCigar(Cigar cigar) {
        return ArtificialSAMUtils.createArtificialRead(Utils.arrayFromArrayWithLength(BASES, cigar.getReadLength()), Utils.arrayFromArrayWithLength(QUALS, cigar.getReadLength()), cigar.toString());
    }

    public static GATKSAMRecord makeReadFromCigar(String cigarString) {
        return makeReadFromCigar(cigarFromString(cigarString));
    }

        /**
        * This function generates every valid permutation of cigar strings with a given length.
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
    public static List<Cigar> generateCigarList(int maximumLength) {
        int numCigarElements = cigarElements.length;
        LinkedList<Cigar> cigarList = new LinkedList<Cigar>();
        byte [] cigarCombination = new byte[maximumLength];

        Utils.fillArrayWithByte(cigarCombination, (byte) 0);               // we start off with all 0's in the combination array.
        int currentIndex = 0;
        while (true) {
            Cigar cigar = createCigarFromCombination(cigarCombination);    // create the cigar
            cigar = combineAdjacentCigarElements(cigar);                   // combine adjacent elements
            if (isCigarValid(cigar)) {                                     // check if it's valid
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

    private static boolean isCigarValid(Cigar cigar) {
        if (cigar.isValid(null, -1) == null) {                                                                          // This should take care of most invalid Cigar Strings (picard's "exhaustive" implementation)

            Stack<CigarElement> cigarElementStack = new Stack<CigarElement>();                                          // Stack to invert cigar string to find ending operator
            CigarOperator startingOp = null;
            CigarOperator endingOp = null;

            // check if it doesn't start with deletions
            boolean readHasStarted = false;                                                                             // search the list of elements for the starting operator
            for (CigarElement cigarElement : cigar.getCigarElements()) {
                if (!readHasStarted) {
                    if (cigarElement.getOperator() != CigarOperator.SOFT_CLIP && cigarElement.getOperator() != CigarOperator.HARD_CLIP) {
                        readHasStarted = true;
                        startingOp = cigarElement.getOperator();
                    }
                }
                cigarElementStack.push(cigarElement);
            }

            while (!cigarElementStack.empty()) {
                CigarElement cigarElement = cigarElementStack.pop();
                if (cigarElement.getOperator() != CigarOperator.SOFT_CLIP && cigarElement.getOperator() != CigarOperator.HARD_CLIP) {
                    endingOp = cigarElement.getOperator();
                    break;
                }
            }

              if (startingOp != CigarOperator.DELETION && endingOp != CigarOperator.DELETION)
                  return true;                                                                                          // we don't accept reads starting or ending in deletions (add any other constraint here)
        }

        return false;
    }

    private static Cigar createCigarFromCombination(byte[] cigarCombination) {
        Cigar cigar = new Cigar();
        for (byte i : cigarCombination) {
            cigar.add(cigarElements[i]);
        }
        return cigar;
    }


    /**
     * Combines equal adjacent elements of a Cigar object
     *
     * @param rawCigar the cigar object
     * @return a combined cigar object
     */
    private static Cigar combineAdjacentCigarElements(Cigar rawCigar) {
        Cigar combinedCigar = new Cigar();
        CigarElement lastElement = null;
        int lastElementLength = 0;
        for (CigarElement cigarElement : rawCigar.getCigarElements()) {
            if (lastElement != null &&
                    ((lastElement.getOperator() == cigarElement.getOperator()) ||
                     (lastElement.getOperator() == CigarOperator.I && cigarElement.getOperator() == CigarOperator.D) ||
                     (lastElement.getOperator() == CigarOperator.D && cigarElement.getOperator() == CigarOperator.I)))
                lastElementLength += cigarElement.getLength();
            else
            {
                if (lastElement != null)
                    combinedCigar.add(new CigarElement(lastElementLength, lastElement.getOperator()));

                lastElement = cigarElement;
                lastElementLength = cigarElement.getLength();
            }
        }
        if (lastElement != null)
            combinedCigar.add(new CigarElement(lastElementLength, lastElement.getOperator()));

        return combinedCigar;
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

    public static Cigar invertCigar (Cigar cigar) {
        Stack<CigarElement> cigarStack = new Stack<CigarElement>();
        for (CigarElement cigarElement : cigar.getCigarElements())
            cigarStack.push(cigarElement);

        Cigar invertedCigar = new Cigar();
        while (!cigarStack.isEmpty())
            invertedCigar.add(cigarStack.pop());

        return invertedCigar;
    }

    /**
     * Checks whether or not the read has any cigar element that is not H or S
     *
     * @param read the read
     * @return true if it has any M, I or D, false otherwise
     */
    public static boolean readHasNonClippedBases(GATKSAMRecord read) {
        for (CigarElement cigarElement : read.getCigar().getCigarElements())
            if (cigarElement.getOperator() != CigarOperator.SOFT_CLIP && cigarElement.getOperator() != CigarOperator.HARD_CLIP)
                return true;
        return false;
    }

    public static Cigar cigarFromString(String cigarString) {
        return TextCigarCodec.getSingleton().decode(cigarString);
    }
}
