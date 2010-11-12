/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.variantrecalibration;

import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.BaseTest;
import org.testng.annotations.Test;
import org.testng.annotations.BeforeTest;
import org.testng.Assert;

import Jama.*; 

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;
import java.util.regex.Pattern;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Feb 26, 2010
 */

public final class VariantGaussianMixtureModelUnitTest extends BaseTest {
    private static int N_VARIANTS = 100;
    VariantDatum[] variantData1 = new VariantDatum[N_VARIANTS];

    @BeforeTest
    public void beforeTest() {
        for ( int i = 0; i < N_VARIANTS; i++ ) {
            variantData1[i].isKnown = i % 2 == 0;      // every other is know
            variantData1[i].qual = (N_VARIANTS - i) * 1.0;
        }

        // first 25 are tv, 25-75 are 50/50, and 75+ are all transitions
        int i = 0;
        for ( ; i < (N_VARIANTS * 0.25); i++ ) { variantData1[i].isTransition = true; }
        for ( ; i < (N_VARIANTS * 0.75); i++ ) { variantData1[i].isTransition = i % 2 == 0; }
        for ( ; i < N_VARIANTS; i++ ) { variantData1[i].isTransition = false; }
    }

    @Test
    public final void testFindTranches1() {
        List<Tranche> tranches = VariantGaussianMixtureModel.findTranches(variantData1, new double[]{0.1, 20}, 2.0);
        Assert.assertEquals( tranches.size(), 2 );

        Tranche t1 = tranches.get(0);
        Assert.assertEquals( t1.fdr, 0.1 );
        Assert.assertEquals( t1.pCut, 26 );
        Assert.assertEquals( t1.numKnown, 37 );
        Assert.assertEquals( t1.numNovel, 37 );

        Tranche t2 = tranches.get(1);
        Assert.assertEquals( t2.fdr, 20 );
        Assert.assertEquals( t2.pCut, 21 );
        Assert.assertEquals( t2.numKnown, 37 );
        Assert.assertEquals( t2.numNovel, 37 );
    }
}