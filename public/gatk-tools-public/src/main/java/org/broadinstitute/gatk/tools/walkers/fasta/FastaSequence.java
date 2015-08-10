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

package org.broadinstitute.gatk.tools.walkers.fasta;

import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

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
        if ( printedHeader ) throw new ReviewedGATKException("Can not set name for FASTA record: header is already printed.");
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