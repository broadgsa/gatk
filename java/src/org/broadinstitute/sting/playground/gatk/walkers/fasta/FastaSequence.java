package org.broadinstitute.sting.playground.gatk.walkers.fasta;

import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.PrintStream;

// fasta sequence holder class

public class FastaSequence {

    private PrintStream out;
    private StringBuffer sb = new StringBuffer();
    private long sequenceCounter = 1;
    private boolean printedHeader = false;

	public FastaSequence(PrintStream out) {
        this.out = out;
    }

    public void append(String s) {
        sb.append(s);
        printFasta(false);
    }

    public void flush() {
        printFasta(true);
        printedHeader = false;
        sequenceCounter++;
    }

    public String getCurrentID() {
        return String.valueOf(sequenceCounter);
    }

    private void printFasta(boolean printAll) {
        if ( sb.length() == 0 || (!printAll && sb.length() < 60) )
            return;
        if ( !printedHeader ) {
            out.println(">" + sequenceCounter);
            printedHeader = true;
        }
        int lines = sb.length() / 60;
        int currentStart = 0;
        for (int i=0; i < lines; i++) {
            out.println(sb.substring(currentStart, currentStart+60));
            currentStart += 60;
        }
        if ( printAll ) {
            out.println(sb.substring(currentStart));
            sb.setLength(0);
        } else {
            sb.delete(0, currentStart);
        }
    }
}