package org.broadinstitute.sting.gatk.walkers.fasta;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.PrintStream;

// fasta sequence holder class

public class FastaSequence {

    private PrintStream out;
    private StringBuffer sb = new StringBuffer();
    private long sequenceCounter = 1;
    private boolean printedHeader = false;
    private String name = null;
    private int lineWidth = -1;
    private boolean noHeader = false;

	public FastaSequence(PrintStream out, int lineWidth, boolean noHeader) {
        this.out = out;
        this.lineWidth = lineWidth;
        this.noHeader = noHeader;
    }

    public void setName(String name) {
        if ( printedHeader ) throw new ReviewedStingException("Can not set name for FASTA record: header is already printed.");
        this.name = name;
    }

    public String getName() {
        if ( name != null ) return name;
        else return getCurrentID();                        
    }
    
    public void append(String s) {
        sb.append(s);
        printFasta(false);
    }

    public void flush() {
        printFasta(true);
        printedHeader = false;
        name = null;
        sequenceCounter++;
    }

    public long getCurrentCount() {
        return sequenceCounter;
    }

    public String getCurrentID() {
        return String.valueOf(sequenceCounter);
    }

    private void printFasta(boolean printAll) {
        if ( sb.length() == 0 || (!printAll && sb.length() < lineWidth) )
            return;
        if ( !printedHeader && !noHeader) {
            if ( name == null ) out.println(">" + sequenceCounter);
            else out.println(">" + name);
            printedHeader = true;
        }
        int lines = sb.length() / lineWidth;
        int currentStart = 0;
        for (int i=0; i < lines; i++) {
            out.println(sb.substring(currentStart, currentStart+lineWidth));
            currentStart += lineWidth;
        }
        if ( printAll ) {
            out.println(sb.substring(currentStart));
            sb.setLength(0);
        } else {
            sb.delete(0, currentStart);
        }
    }
}