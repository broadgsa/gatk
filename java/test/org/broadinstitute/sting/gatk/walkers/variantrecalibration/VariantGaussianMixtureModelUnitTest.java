/*
 * Copyright (c) 2010, The Broad Institute
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
import org.broadinstitute.sting.utils.exceptions.StingException;
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

    private final File QUAL_DATA = new File(testDir + "tranches.raw.dat");
    private final double[] FDRS = new double[]{0.1, 1, 10, 100};
    private final double TARGET_TITV = 2.8;
    private final File EXPECTED_TRANCHES_NEW = new File(testDir + "tranches.6.txt");
    private final File EXPECTED_TRANCHES_OLD = new File(testDir + "tranches.4.txt");

    private List<VariantDatum> readData() {
        List<VariantDatum> vd = new ArrayList<VariantDatum>();
        try {
            for ( String line : new XReadLines(QUAL_DATA, true) ) {
                String[] parts = line.split(" ");
                if ( ! parts[0].equals("QUAL") ) {
                    VariantDatum datum = new VariantDatum();
                    datum.qual = Double.valueOf(parts[0]);
                    datum.isTransition = parts[1].equals("1");
                    datum.isKnown = parts[2].equals("1");
                    vd.add(datum);
                }
            }
        } catch (FileNotFoundException e) {
            throw new StingException("foo", e);
        }

        return vd;
    }

    @Test(expectedExceptions = {UserException.MalformedFile.class})
    public final void readBadFormat() {
        Tranche.readTraches(QUAL_DATA);
    }

    @Test
    public final void readNewFormat() {
        read(EXPECTED_TRANCHES_NEW);
    }

    @Test
    public final void readOldFormat() {
        read(EXPECTED_TRANCHES_OLD);
    }

    public final List<Tranche> read(File f) {
        return Tranche.readTraches(f);
    }

    @Test
    public final void testNewAndOldAretheSame() {
        List<Tranche> newFormat = read(EXPECTED_TRANCHES_NEW);
        List<Tranche> oldFormat = read(EXPECTED_TRANCHES_OLD);
        assertTranchesAreTheSame(newFormat, oldFormat, false, true);
    }

    private static void assertTranchesAreTheSame(List<Tranche> newFormat, List<Tranche> oldFormat, boolean completeP, boolean includeName) {
        Assert.assertEquals(oldFormat.size(), newFormat.size());
        for ( int i = 0; i < newFormat.size(); i++ ) {
            Tranche n = newFormat.get(i);
            Tranche o = oldFormat.get(i);
            Assert.assertEquals(n.fdr, o.fdr, 1e-3);
            Assert.assertEquals(n.numNovel, o.numNovel);
            Assert.assertEquals(n.novelTiTv, o.novelTiTv, 1e-3);
            if ( includeName )
                Assert.assertEquals(n.name, o.name);
            if ( completeP ) {
                Assert.assertEquals(n.numKnown, o.numKnown);
                Assert.assertEquals(n.knownTiTv, o.knownTiTv, 1e-3);
                Assert.assertEquals(n.targetTiTv, o.targetTiTv, 1e-3);
            }
        }
    }

    @Test
    public final void testFindTranches1() {
        List<VariantDatum> vd = readData();
        List<Tranche> tranches = VariantGaussianMixtureModel.findTranches(vd.toArray(new VariantDatum[0]), FDRS, TARGET_TITV);
        System.out.printf(Tranche.tranchesString(tranches));
        assertTranchesAreTheSame(read(EXPECTED_TRANCHES_NEW), tranches, true, false);
    }

//    @Test(expectedExceptions = {UserException.class})
//    public final void testBadFDR() {
//        List<VariantDatum> vd = readData();
//        VariantGaussianMixtureModel.findTranches(vd.toArray(new VariantDatum[0]), new double[]{-1}, TARGET_TITV);
//    }
//
//    @Test(expectedExceptions = {UserException.class})
//    public final void testBadTargetTiTv() {
//        List<VariantDatum> vd = readData();
//        VariantGaussianMixtureModel.findTranches(vd.toArray(new VariantDatum[0]), FDRS, 0.1);
//    }
}
