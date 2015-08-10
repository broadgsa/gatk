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

package org.broadinstitute.gatk.engine.io;

import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.ReadUtils;

import java.io.FileNotFoundException;
import java.io.PrintStream;

/**
 * User: carneiro
 * Date: 1/27/13
 * Time: 12:54 AM
 */
public class FastqFileWriter {
    private PrintStream output;

    public FastqFileWriter(String filename) {
        try {
            this.output = new PrintStream(filename);
        } catch (FileNotFoundException e) {
            throw new ReviewedGATKException("Can't open file " + filename);
        }
    }

    public void addAlignment(GATKSAMRecord read) {
        output.println("@" + read.getReadName());

        if (read.getReadNegativeStrandFlag()) {
            output.println(ReadUtils.getBasesReverseComplement(read));
            output.println("+");
            output.println(ReadUtils.convertReadQualToString(invertQuals(read.getBaseQualities())));
        } else {
            output.println(ReadUtils.convertReadBasesToString(read));
            output.println("+");
            output.println(ReadUtils.convertReadQualToString(read));
        }
    }

    public void close() {
        this.output.close();
    }

    private byte[] invertQuals (byte[] quals) {
        final int l = quals.length;
        byte[] invertedQuals = new byte[l];
        for (int i=0; i<l; i++) {
            invertedQuals[l-1-i] = quals[i];
        }
        return invertedQuals;
    }
}
