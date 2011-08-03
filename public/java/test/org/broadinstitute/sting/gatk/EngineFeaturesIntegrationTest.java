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

package org.broadinstitute.sting.gatk;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.List;

/**
 *
 */
public class EngineFeaturesIntegrationTest extends WalkerTest {
    private void testBadRODBindingInput(String type, String name) {
        WalkerTestSpec spec = new WalkerTestSpec("-T SelectVariants -L 1:1 --variants:" + type + " "
                + b37dbSNP132 + " -R " + b37KGReference + " -o %s",
                1, UserException.class);
        executeTest(name, spec);
    }


    @Test() private void testBadRODBindingInputType1() {
        testBadRODBindingInput("beagle", "BEAGLE input to VCF expecting walker");
    }

    @Test() private void testBadRODBindingInputType2() {
        testBadRODBindingInput("vcf3", "VCF3 input to VCF expecting walker");
    }

    @Test() private void testBadRODBindingInputType3() {
        testBadRODBindingInput("bed", "Bed input to VCF expecting walker");
    }

    @Test() private void testBadRODBindingInputTypeUnknownType() {
        testBadRODBindingInput("bedXXX", "Unknown input to VCF expecting walker");
    }
}

//class TestRodBindings extends RodWalker<Integer, Integer> {
//    @Input(fullName="req", required=true)
//    public RodBinding<Feature> required;
//
//    @Input(fullName="optional", required=false)
//    public RodBinding<Feature> optional = RodBinding.makeUnbound(Feature.class);
//
//    @Input(fullName="rodList", shortName="RL", doc="A list of ROD types that we will convert to a table", required=true)
//    public List<RodBinding<Feature>> variantsList;
//
//    public void initialize() {
//        // bound values
//        Assert.assertEquals(required.isBound(), true);
//
//
//        System.exit(0);
//    }
//
//    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) { return 0; }
//    public Integer reduceInit() { return 0; }
//    public Integer reduce(Integer counter, Integer sum) { return counter + sum; }
//}