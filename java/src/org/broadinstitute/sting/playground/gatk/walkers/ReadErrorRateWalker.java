package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.BaseUtils;
import net.sf.samtools.SAMRecord;

import java.util.Random;

/**
 * ReadErrorRateWalker assesses the error rate per read position ('cycle') by comparing the
 * read to its home on the reference and noting the mismatch rate.  It ignores reads with
 * indels in them, treats high and low-quality reference bases the same, and does not count
 * ambiguous bases as mismatches.  It's also thread-safe, so you can process a slew of reads
 * in short order.
 *
 * @author Kiran Garimella
 */
public class ReadErrorRateWalker extends ReadWalker<boolean[], int[]> {
    @Argument(fullName="printVisualHits",   shortName="v",  doc="print visual hits",    required=false) public boolean printVisualHits = false;
    @Argument(fullName="useNextBestBase",   shortName="nb", doc="use next best base",   required=false) public boolean useNextBestBase = false;
    @Argument(fullName="useNextRandomBase", shortName="nr", doc="use next random base", required=false) public boolean useNextRandomBase = false;

    /**
     * Ignore reads with indels or clipping
     *
     * @param read     the read to assess
     * @return true if the read can be processed, false if it should be ignored
     */
    public boolean filter(char[] ref, SAMRecord read) {
        return (read.getCigar().numCigarElements() == 1 && read.getReadLength() <= ref.length);
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
    public boolean[] map(char[] ref, SAMRecord read) {
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
                byte compBase = convertIUPACBaseToSimpleBase((byte)ref[cycle]);

                System.out.print((char) compBase);
            }
            System.out.println("\n");
        }

        for (int cycle = 0; cycle < bases.length; cycle++) {
            byte compBase = convertIUPACBaseToSimpleBase((byte)ref[cycle]);

            if (compBase != '.') {
                if (useNextBestBase || useNextRandomBase) {
                    byte nextBestBase;
                    if (useNextBestBase) {
                        nextBestBase = (byte) BaseUtils.baseIndexToSimpleBase(QualityUtils.compressedQualityToBaseIndex(sq[cycle]));
                    } else {
                        nextBestBase = bases[cycle];
                        Random generator = new Random();
                        while (nextBestBase == bases[cycle]) {
                            nextBestBase = (byte) BaseUtils.baseIndexToSimpleBase(generator.nextInt(4));
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
    public int[] reduceInit() {
        return null;
    }

    /**
     * Summarize the error rate data.
     *
     * @param value  the read mismatch array
     * @param sum    the summed mismatch array
     * @return the summed mismatch array with the new read mismatch array added
     */
    public int[] reduce(boolean[] value, int[] sum) {
        if (sum == null) {
            sum = new int[value.length];
        }

        for (int cycle = 0; cycle < value.length; cycle++) {
            sum[cycle] += (value[cycle] ? 1 : 0);
        }

        return sum;
    }

    /**
     * We've processed all the reads.  Emit the final, normalized error rate data.
     *
     * @param sum  the summed mismatch array
     */
    public void onTraversalDone(int[] sum) {
        for (int cycle = 0; cycle < sum.length - 1; cycle++) {
            double errorrate = ((double) sum[cycle])/((double) sum[sum.length - 1]);
            out.println("[ERROR_RATE] " + cycle + " " + errorrate);
        }
    }
}
