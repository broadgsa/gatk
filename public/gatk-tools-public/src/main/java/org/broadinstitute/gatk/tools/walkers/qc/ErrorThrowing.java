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

package org.broadinstitute.gatk.tools.walkers.qc;

import org.broadinstitute.gatk.engine.walkers.FailMethod;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.NanoSchedulable;
import org.broadinstitute.gatk.engine.walkers.RefWalker;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;

/**
 * A walker that simply throws errors.  Allows us to test that the engine is behaving as expected with error handling
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_TOY, extraDocs = {CommandLineGATK.class} )
public class ErrorThrowing extends RefWalker<Integer,Integer> implements TreeReducible<Integer>, NanoSchedulable {
    @Input(fullName="exception", shortName = "E", doc="Java class of exception to throw", required=true)
    public String exceptionToThrow;

    @Argument(fullName = "failMethod", shortName = "fail", doc = "Determines which method to fail in", required = false)
    public FailMethod failMethod = FailMethod.MAP;

    //
    // Template code to allow us to build the walker, doesn't actually do anything
    //
    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( ref == null ) // only throw exception when we are in proper map, not special map(null) call
            return null;

        if ( failMethod == FailMethod.MAP )
            FailMethod.fail(exceptionToThrow);

        return 0;
    }

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        if ( value != null && failMethod == FailMethod.REDUCE )
            FailMethod.fail(exceptionToThrow);
        return sum;
    }

    public Integer treeReduce(final Integer lhs, final Integer rhs) {
        if ( failMethod == FailMethod.TREE_REDUCE )
            FailMethod.fail(exceptionToThrow);
        return rhs;
    }
}
