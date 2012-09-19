/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.qc;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.NanoSchedulable;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;

/**
 * a walker that simply throws errors.  Allows us to test that the engine is behaving as expected with error handling
 */
@Hidden
@DocumentedGATKFeature( groupName = "Quality Control and Simple Analysis Tools", extraDocs = {CommandLineGATK.class} )
public class ErrorThrowing extends RefWalker<Integer,Integer> implements TreeReducible<Integer>, NanoSchedulable {
    @Input(fullName="exception", shortName = "E", doc="Java class of exception to throw", required=true)
    public String exceptionToThrow;

    @Argument(fullName = "failMethod", shortName = "fail", doc = "Determines which method to fail in", required = false)
    public FailMethod failMethod = FailMethod.MAP;

    public enum FailMethod {
          MAP,
          REDUCE,
          TREE_REDUCE
    }

    //
    // Template code to allow us to build the walker, doesn't actually do anything
    //
    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( ref == null ) // only throw exception when we are in proper map, not special map(null) call
            return null;

        if ( failMethod == FailMethod.MAP )
            fail();

        return 0;
    }

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        if ( value != null && failMethod == FailMethod.REDUCE )
            fail();
        return sum;
    }

    public Integer treeReduce(final Integer lhs, final Integer rhs) {
        if ( failMethod == FailMethod.TREE_REDUCE )
            fail();
        return rhs;
    }

    private void fail() {
        if ( exceptionToThrow.equals("UserException") ) {
            throw new UserException("UserException");
        } else if ( exceptionToThrow.equals("NullPointerException") ) {
            throw new NullPointerException();
        } else if ( exceptionToThrow.equals("ReviewedStingException") ) {
            throw new ReviewedStingException("ReviewedStingException");
        } else if ( exceptionToThrow.equals("SamError1") ) {
            throw new RuntimeException(CommandLineGATK.PICARD_TEXT_SAM_FILE_ERROR_1);
        } else if ( exceptionToThrow.equals("SamError2") ) {
            throw new RuntimeException(CommandLineGATK.PICARD_TEXT_SAM_FILE_ERROR_2);
        } else if ( exceptionToThrow.equals("NoSpace1") ) {
            throw new net.sf.samtools.util.RuntimeIOException(new java.io.IOException("No space left on device java.io.FileOutputStream.writeBytes(Native Method)"));
        } else if ( exceptionToThrow.equals("NoSpace2") ) {
            throw new net.sf.samtools.SAMException("Exception writing BAM index file", new java.io.IOException("No space left on device java.io.FileOutputStream.writeBytes(Native Method)"));
        } else {
            throw new UserException.BadArgumentValue("exception", "exception isn't a recognized value " + exceptionToThrow);
        }
    }
}
