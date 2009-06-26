package org.broadinstitute.sting.utils.genotype;

import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.assertEquals;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.genotype.LikelihoodObject;
import static junit.framework.Assert.assertTrue;


/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @author aaron
 *         <p/>
 *         Class LikelihoodObjectTest
 *         <p/>
 *         Tests the Likelihood object.
 */
public class LikelihoodObjectTest extends BaseTest {

    private LikelihoodObject mLO = null;

    @Before
    public void before() {
        mLO = new LikelihoodObject();
    }

    @Test
    public void testBlankConstruction() {
        mLO = new LikelihoodObject();
        assertTrue(mLO.likelihood.size() == LikelihoodObject.GENOTYPE.values().length);
    }

    @Test
    public void testConstructionFromArray() {
        double[] ray = new double[10];
        for (int x = 0; x < 10; x++) {
            ray[x] = ( x * 25 );
        }
        mLO = new LikelihoodObject(ray);
        assertTrue(mLO.likelihood.size() == LikelihoodObject.GENOTYPE.values().length);

        int index = 0;
        for (LikelihoodObject.GENOTYPE t : LikelihoodObject.GENOTYPE.values()) {
            assertTrue(ray[index] == mLO.likelihood.get(t));
            ++index;
        }
    }

    @Test
    public void testByteArrayReturn() {
        double[] ray = new double[10];
        for (int x = 0; x < 10; x++) {
            ray[x] = ( x * 25.0 );
        }
        mLO = new LikelihoodObject(ray);
        assertTrue(mLO.likelihood.size() == LikelihoodObject.GENOTYPE.values().length);

        int index = 0;
        short[] ret = mLO.toByteArray();
        for (index = 0; index < ret.length; index++) {
            assertTrue(ray[index] == ret[index]);
        }
    }

    @Test
    public void testDefaultArrayValues() {
        mLO = new LikelihoodObject();
        short[] ret = mLO.toByteArray();
        for (int index = 0; index < ret.length; index++) {
            assertTrue(ret[index] == 255);
        }
    }

    @Test
    public void testGetMinimum() {
        double[] ray = new double[10];
        for (int x = 0; x < 10; x++) {
            ray[x] = ( 240.0 );
            ray[x] = ( 240.0 );
        }
        ray [5] = 0;
        mLO = new LikelihoodObject(ray);
        assertTrue(mLO.likelihood.size() == LikelihoodObject.GENOTYPE.values().length);
        short smallest = (short)mLO.getBestLikelihood();
        assertTrue(smallest == 0);
        int index = 0;
        short[] ret = mLO.toByteArray();
        for (index = 0; index < ret.length; index++) {
            assertTrue(smallest <= ret[index]);
        }
    }


    @Test
    public void testSetLikelihood() {
        mLO = new LikelihoodObject();
        for (LikelihoodObject.GENOTYPE t : LikelihoodObject.GENOTYPE.values()) {
            mLO.setLikelihood(t,128);
        }
        assertTrue(mLO.likelihood.size() == LikelihoodObject.GENOTYPE.values().length);

        int index = 0;
        short[] ret = mLO.toByteArray();
        for (index = 0; index < ret.length; index++) {
            assertTrue(ret[index] == 128);
        }
    }


}
