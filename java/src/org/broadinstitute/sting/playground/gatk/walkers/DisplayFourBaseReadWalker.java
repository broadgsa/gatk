package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.QualityUtils;
import net.sf.samtools.SAMRecord;

public class DisplayFourBaseReadWalker extends ReadWalker<Integer, Integer> {
    public Integer map(LocusContext context, SAMRecord read) {
        String bases = read.getReadString();

        boolean displayed = false;

        byte[] sq = (byte[]) read.getAttribute("SQ");

        if (read.getReadName().equalsIgnoreCase("30JJE.5.24197751")) {
            System.out.println(read.format());
        }

        for (int i = 0; i < sq.length; i++) {
            int baseIndex = QualityUtils.compressedQualityToBaseIndex(sq[i]);
            char base = '.';
            switch (baseIndex) {
                case 0: base = 'A'; break;
                case 1: base = 'C'; break;
                case 2: base = 'G'; break;
                case 3: base = 'T'; break;
            }

            if (base == bases.charAt(i)) {
                if (!displayed) {
                    System.out.println(bases);
                    displayed = true;
                }
                
                System.out.print(base);
            }
            else {
                if (displayed) {
                    System.out.print(" ");
                }
            }
        }
        if (displayed) {
            System.out.print("\n");
        }

        return 0;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return 0;
    }
}
