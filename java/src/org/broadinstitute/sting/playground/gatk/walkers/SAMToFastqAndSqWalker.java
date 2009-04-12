package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import net.sf.samtools.SAMRecord;

import java.io.*;

public class SAMToFastqAndSqWalker extends ReadWalker<Integer, Integer> {
    @Argument(fullName="FASTQFILE", doc="Output path for Fastq file")
    public File FASTQFILE;

    @Argument(fullName="SQFILE", shortName="SQ", doc="Output path for secondary quality map (readName => SAM SQ field)", required=false)
    public File SQFILE;

    private PrintWriter fastqbw;
    private PrintWriter sqbw;

    private boolean ready = false;
    
    public Integer map(LocusContext context, SAMRecord read) {
        if (!ready) {
            try {
                fastqbw = new PrintWriter(FASTQFILE);

                if (SQFILE != null) {
                    sqbw = new PrintWriter(SQFILE);
                }

                ready = true;
            } catch (IOException e) {
                err.println("Unable to open output files.");
                System.exit(1);
            }
        }

        String bases = read.getReadString();
        String quals = read.getBaseQualityString();
        byte[] sqs   = (byte[]) read.getAttribute("SQ");

        fastqbw.println("@" + read.getReadName());
        fastqbw.println(bases);

        fastqbw.println("+" + read.getReadName());
        fastqbw.println(quals);

        if (sqbw != null && sqs != null) {
            sqbw.print(read.getReadName() + "\t" + "SQ:H:");

            for (byte sq : sqs) {
                sqbw.printf("%02X", sq);
            }
            sqbw.println();
        }

        return null;
    }

    public Integer reduceInit() {
        return null;
    }

    public Integer reduce(Integer value, Integer sum) {
        return null;
    }

    public void onTraversalDone(Integer sum) {
        fastqbw.close();

        if (sqbw != null) {
            sqbw.close();
        }
    }
}
