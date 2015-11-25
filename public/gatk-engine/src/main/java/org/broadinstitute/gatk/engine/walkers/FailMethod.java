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

package org.broadinstitute.gatk.engine.walkers;

import htsjdk.samtools.SAMException;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;

public enum FailMethod {
      MAP,
      REDUCE,
      TREE_REDUCE;

    /**
     * Used by QC walkers to test that engine throws appropriate errors.
     * Split from the walker in ErrorThrowing.java.
     * @param exceptionToThrow Exception type to throw.
     */
    public static void fail(final String exceptionToThrow) {
        switch (exceptionToThrow) {
            case "UserException":
                throw new UserException("UserException");
            case "NullPointerException":
                throw new NullPointerException();
            case "ReviewedGATKException":
                throw new ReviewedGATKException("ReviewedGATKException");
            case "SamError1":
                throw new RuntimeException(CommandLineGATK.PICARD_TEXT_SAM_FILE_ERROR_1);
            case "SamError2":
                throw new RuntimeException(CommandLineGATK.PICARD_TEXT_SAM_FILE_ERROR_2);
            case "NoSpace1":
                throw new htsjdk.samtools.util.RuntimeIOException(new java.io.IOException("No space left on device java.io.FileOutputStream.writeBytes(Native Method)"));
            case "NoSpace2":
                throw new SAMException("Exception writing BAM index file", new java.io.IOException("No space left on device java.io.FileOutputStream.writeBytes(Native Method)"));
            default:
                throw new UserException.BadArgumentValue("exception", "exception isn't a recognized value " + exceptionToThrow);
        }
    }
}
