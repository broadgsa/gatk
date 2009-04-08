package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.QualityUtils;
import net.sf.samtools.SAMRecord;
import edu.mit.broad.picard.reference.ReferenceSequence;

/**
 * ReadErrorRateWalker assesses the error rate per read position ('cycle') by comparing
 */
public class ReadErrorRateWalker extends ReadWalker<boolean[], int[]> {
    @Argument(fullName="useNextBestBase",required=false,defaultValue="false")
    public boolean useNextBestBase;

    public boolean filter(LocusContext context, SAMRecord read) {
        return (read.getCigar().numCigarElements() == 1);
    }

    public boolean[] map(LocusContext context, SAMRecord read) {
        boolean[] errorsPerCycle = new boolean[read.getReadLength() + 1];

        byte[] bases  = read.getReadBases();
        byte[] contig = context.getReferenceContig().getBases();
        byte[] sq     = (byte[]) read.getAttribute("SQ");

        int totalMismatches = 0;

        for (int cycle = 0, offset = (int) context.getPosition(); cycle < bases.length; cycle++, offset++) {
            byte compBase;

            switch (contig[offset]) {
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

            if (compBase != '.') {
                if (useNextBestBase) {
                    int nextBestBaseIndex = QualityUtils.compressedQualityToBaseIndex(sq[cycle]);
                    byte nextBestBase;
                    switch (nextBestBaseIndex) {
                        case 0: nextBestBase = 'A'; break;
                        case 1: nextBestBase = 'C'; break;
                        case 2: nextBestBase = 'G'; break;
                        case 3: nextBestBase = 'T'; break;
                        default: nextBestBase = '.'; break;
                    }

                    if (nextBestBase != '.') {
                        errorsPerCycle[cycle] = !(bases[cycle] == compBase || bases[cycle] == nextBestBase);
                        totalMismatches = !(bases[cycle] == compBase || bases[cycle] == nextBestBase) ? 1 : 0;
                    }
                } else {
                    errorsPerCycle[cycle] = (bases[cycle] != compBase);
                    totalMismatches += (bases[cycle] != compBase) ? 1 : 0;
                }
            }
        }

        if (totalMismatches > 4) {
            for (int cycle = 0; cycle < bases.length; cycle++) { System.out.print((char) bases[cycle]); } System.out.print("\n");
            for (int cycle = 0, offset = (int) context.getPosition(); cycle < bases.length; cycle++, offset++) { System.out.print((char) contig[offset]); } System.out.print("\n");
            System.out.println(totalMismatches + "\n");
        }

        // We encode that we saw a read in the last position of the array.
        // That way we know what to normalize by, and we get thread safety!
        errorsPerCycle[errorsPerCycle.length - 1] = true;

        return errorsPerCycle;
    }

    public int[] reduceInit() {
        return null;
    }

    public int[] reduce(boolean[] value, int[] sum) {
        if (sum == null) {
            sum = new int[value.length];
        }

        for (int cycle = 0; cycle < value.length; cycle++) {
            sum[cycle] += (value[cycle] ? 1 : 0);
        }

        return sum;
    }

    public void onTraversalDone(int[] sum) {
        for (int cycle = 0; cycle < sum.length - 1; cycle++) {
            double errorrate = ((double) sum[cycle])/((double) sum[sum.length - 1]);
            System.out.println("[ERROR_RATE] " + cycle + " " + errorrate);
        }
    }

}
