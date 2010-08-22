/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.BaseUtils;
import net.sf.samtools.SAMRecord;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.io.PrintStream;

/**
 * ReadErrorRateWalker assesses the error rate per read position ('cycle') by comparing the
 * read to its home on the reference and noting the mismatch rate.  It ignores reads with
 * indels in them, treats high and low-quality reference bases the same, and does not count
 * ambiguous bases as mismatches.  It's also thread-safe, so you can process a slew of reads
 * in short order.
 *
 * @author Kiran Garimella
 */
public class ReadErrorRateWalker extends ReadWalker<boolean[], ReadErrorRateCollection> implements TreeReducible<ReadErrorRateCollection> {
    @Output PrintStream out;
    @Argument(fullName="printVisualHits",   shortName="v",  doc="print visual hits",    required=false) public boolean printVisualHits = false;
    @Argument(fullName="useNextBestBase",   shortName="nb", doc="use next best base",   required=false) public boolean useNextBestBase = false;
    @Argument(fullName="useNonNextBestBase",shortName="nnb",doc="use nonnext best base",required=false) public boolean useNonNextBestBase = false;
    @Argument(fullName="useNextRandomBase", shortName="nr", doc="use next random base", required=false) public boolean useNextRandomBase = false;

    /**
     * Ignore reads with indels or clipping
     *
     * @param read     the read to assess
     * @return true if the read can be processed, false if it should be ignored
     */
    public boolean filter(ReferenceContext ref, SAMRecord read) {
        return (read.getCigar().numCigarElements() == 1 && read.getReadLength() <= ref.getBases().length && (!useNonNextBestBase || read.getAttribute("SQ") != null));
    }

    /**
     * For each read, return a boolean array indicating the locations of the mismatch.
     * Length of the array is one element longer than the read length.  The last element
     * of this array is always "true" so that we can figure out how many reads we
     * processed in a thread-safe manner.
     * 
     * @param read     the read to assess
     * @return An array of length (read_length + 1) indicating where the mismatches occur.
     *         Last element is for internal use so the reduce() function can figure out how
     *         many reads we processed.
     */
    public boolean[] map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        boolean[] errorsPerCycle = new boolean[read.getReadLength() + 1];

        byte[] bases  = read.getReadBases();
        byte[] sq     = (byte[]) read.getAttribute("SQ");

        if (printVisualHits) {
            System.out.println(read.getReadName());
            for (int cycle = 0; cycle < bases.length; cycle++) {
                System.out.print((char) bases[cycle]);
            }
            System.out.println();

            for (int cycle = 0; cycle < bases.length; cycle++) {
                byte compBase = convertIUPACBaseToSimpleBase(ref.getBases()[cycle]);

                System.out.print((char) compBase);
            }
            System.out.println("\n");
        }

        for (int cycle = 0; cycle < bases.length; cycle++) {
            byte compBase = convertIUPACBaseToSimpleBase(ref.getBases()[cycle]);

            if (compBase != '.') {
                if (useNextBestBase || useNextRandomBase || useNonNextBestBase) {
                    byte nextBestBase;
                    if (useNextBestBase) {
                        nextBestBase = BaseUtils.baseIndexToSimpleBase(QualityUtils.compressedQualityToBaseIndex(sq[cycle]));
                    } else if (useNonNextBestBase) {
                        nextBestBase = bases[cycle];
                        Random generator = new Random();
                        while (nextBestBase == bases[cycle] || nextBestBase == BaseUtils.baseIndexToSimpleBase(QualityUtils.compressedQualityToBaseIndex(sq[cycle]))) {
                            nextBestBase = BaseUtils.baseIndexToSimpleBase(generator.nextInt(4));
                        }
                    } else {
                        nextBestBase = bases[cycle];
                        Random generator = new Random();
                        while (nextBestBase == bases[cycle]) {
                            nextBestBase = BaseUtils.baseIndexToSimpleBase(generator.nextInt(4));
                        }
                    }

                    if (nextBestBase != '.') {
                        if (read.getReadNegativeStrandFlag()) {
                            errorsPerCycle[bases.length - cycle - 1] = !(bases[cycle] == compBase || nextBestBase == compBase);
                        } else {
                            errorsPerCycle[cycle] = !(bases[cycle] == compBase || nextBestBase == compBase);
                        }
                    }
                } else {
                    if (read.getReadNegativeStrandFlag()) {
                        errorsPerCycle[bases.length - cycle - 1] = !(bases[cycle] == compBase);
                    } else {
                        errorsPerCycle[cycle] = !(bases[cycle] == compBase);
                    }
                }
            }
        }

        // We encode that we saw a read in the last position of the array.
        // That way we know what to normalize by, and we get thread safety!
        errorsPerCycle[errorsPerCycle.length - 1] = true;

        return errorsPerCycle;
    }

    private byte convertIUPACBaseToSimpleBase(byte iupacBase) {
        char compBase;

        switch (iupacBase) {
            case 'A':
            case 'a': compBase = 'A'; break;
            case 'C':
            case 'c': compBase = 'C'; break;
            case 'G':
            case 'g': compBase = 'G'; break;
            case 'T':
            case 't': compBase = 'T'; break;
            default:  compBase = '.'; break;
        }

        return (byte) compBase;
    }

    /**
     * We don't initialize the array here because we need to know how long the read is first.
     *
     * @return null
     */
    public ReadErrorRateCollection reduceInit() {
        return new ReadErrorRateCollection();
    }

    /**
     * Summarize the error rate data.
     *
     * @param value  the read mismatch array
     * @param collection    the summed mismatch array
     * @return the summed mismatch array with the new read mismatch array added
     */
    public ReadErrorRateCollection reduce(boolean[] value, ReadErrorRateCollection collection) {

        collection.update(value);

        return collection;
    }

    /**
     * For multithreading - take two read error rate collections and put them together
     * @param left   one collection
     * @param right   another collection
     * @return left updated with the counts from right
     */
    public ReadErrorRateCollection treeReduce(ReadErrorRateCollection left, ReadErrorRateCollection right) {
        left.merge(right);
        return left;
    }

    /**
     * We've processed all the reads.  Emit the final, normalized error rate data.
     *
     * @param collection  the summed mismatch arrays
     */
    public void onTraversalDone(ReadErrorRateCollection collection) {

        out.print(collection.toString());
    }
}

class ReadErrorRateCollection {
    private HashMap<Integer,int[]> readsByReadLength;

    public ReadErrorRateCollection() {
        readsByReadLength = new HashMap<Integer, int[]>();
    }

    public void update(boolean[] mismatchArray) {
        if ( ! readsByReadLength.containsKey(mismatchArray.length) ) {
            readsByReadLength.put(mismatchArray.length, zeroArray(mismatchArray.length));
        }

        updateErrorCounts(readsByReadLength.get(mismatchArray.length), mismatchArray);
    }

    public String toString() {
        StringBuilder builder = new StringBuilder();
        for ( int length : readsByReadLength.keySet() ) {
            for ( int cycle = 0; cycle < length-1; cycle++) {
                int[] counts = readsByReadLength.get(length);
                builder.append(length);
                builder.append("\t");
                builder.append(cycle);
                builder.append("\t");
                builder.append( ( ( double ) counts[cycle] / ( (double) counts[length-1])));
                builder.append("\n");
            }
        }
        return builder.toString();
    }

    public void merge(ReadErrorRateCollection other) {
        for ( Map.Entry<Integer,int[]> errorCounts : other.readsByReadLength.entrySet() ) {
            if ( this.readsByReadLength.keySet().contains(errorCounts.getKey()) ) {
                mergeCounts(readsByReadLength.get(errorCounts.getKey()),errorCounts.getValue());
            } else {
                readsByReadLength.put(errorCounts.getKey(),errorCounts.getValue());
            }
        }
    }

    private static int[] zeroArray( int length ) {
        int[] array = new int[length];
        for ( int ii = 0; ii < length; ii ++ ) {
            array[ii] = 0;
        }

        return array;
    }

    private static void mergeCounts ( int[] addToMe, int[] dontTouchMe ) {
        for ( int index = 0; index < addToMe.length; index ++ ) {
            addToMe[index] += dontTouchMe[index];
        }
    }

    public static void updateErrorCounts(int[] sum, boolean[] value) {

        for (int cycle = 0; cycle < value.length; cycle++) {
            sum[cycle] += (value[cycle] ? 1 : 0);
        }

    }
}