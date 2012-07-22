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

import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

/**
 * a walker that simply throws errors.  Allows us to test that the engine is behaving as expected with error handling
 */
public class ErrorThrowing extends RodWalker<Integer,Integer> implements TreeReducible<Integer> {
    @Input(fullName="exception", shortName = "E", doc="Java class of exception to throw", required=true)
    public String exceptionToThrow;

    //
    // Template code to allow us to build the walker, doesn't actually do anything
    //
    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( exceptionToThrow.equals("UserException") ) {
            throw new UserException("UserException");
        } else if ( exceptionToThrow.equals("NullPointerException") ) {
            throw new NullPointerException();
        } else if ( exceptionToThrow.equals("ReviewedStingException") ) {
            throw new ReviewedStingException("ReviewedStingException");
        } else {
            throw new UserException.BadArgumentValue("exception", "exception isn't a recognized value " + exceptionToThrow);
        }
    }

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    public Integer treeReduce(final Integer lhs, final Integer rhs) {
        return lhs + rhs;
    }
}
