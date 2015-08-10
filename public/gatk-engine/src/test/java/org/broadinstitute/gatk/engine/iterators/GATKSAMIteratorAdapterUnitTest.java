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

package org.broadinstitute.gatk.engine.iterators;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import org.broadinstitute.gatk.utils.BaseTest;
import static org.testng.Assert.assertEquals;

import org.broadinstitute.gatk.utils.iterators.GATKSAMIterator;
import org.broadinstitute.gatk.utils.iterators.GATKSAMIteratorAdapter;
import org.testng.annotations.Test;

import java.util.Iterator;

/**
 *
 * User: aaron
 * Date: May 13, 2009
 * Time: 6:58:21 PM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 * @date May 13, 2009
 * <p/>
 * Class GATKSAMIteratorTest
 * <p/>
 * Tests the GATKSAMIteratorAdapter class.
 */
public class GATKSAMIteratorAdapterUnitTest extends BaseTest {

    class MyTestIterator implements Iterator<SAMRecord> {

        public int count = 0;

        public MyTestIterator() {
            count = 0;
        }

        public boolean hasNext() {
            if (count < 100) {
                ++count;
                return true;
            } else {
                return false;
            }
        }

        public SAMRecord next() {
            return null;
        }

        public void remove() {
            throw new UnsupportedOperationException("Unsupported");
        }
    }

    class MyTestCloseableIterator implements CloseableIterator<SAMRecord> {
        public int count = 0;

        public MyTestCloseableIterator() {
            count = 0;
        }

        public boolean hasNext() {
            if (count < 100) {
                ++count;
                return true;
            } else {
                return false;
            }
        }

        public SAMRecord next() {
            return null;
        }

        public void remove() {
            throw new UnsupportedOperationException("Unsupported");
        }

        public void close() {
            count = -1;
        }
    }


    @Test
    public void testNormalIterator() {
        final int COUNT = 100;
        MyTestIterator it = new MyTestIterator();

        GATKSAMIterator samIt = GATKSAMIteratorAdapter.adapt(it);
        int countCheck = 0;
        while (samIt.hasNext()) {
            samIt.next();
            ++countCheck;
            //logger.warn("cnt = " + countCheck);
        }

        assertEquals(countCheck, COUNT);

        assertEquals(countCheck, COUNT);
    }

    @Test
    public void testCloseableIterator() {
        final int COUNT = 100;

        MyTestCloseableIterator it = new MyTestCloseableIterator();

        GATKSAMIterator samIt = GATKSAMIteratorAdapter.adapt(it);

        int countCheck = 0;
        while (samIt.hasNext()) {
            samIt.next();
            ++countCheck;
        }

        assertEquals(countCheck, COUNT);
    }

    @Test
    public void testCloseOnCloseableIterator() {
        final int COUNT = 100;

        MyTestCloseableIterator it = new MyTestCloseableIterator();
        
        GATKSAMIterator samIt = GATKSAMIteratorAdapter.adapt(it);


        int countCheck = 0;
        while (samIt.hasNext()) {
            samIt.next();
            ++countCheck;
        }

        assertEquals(countCheck, COUNT);

        // check to see that the count get's set to -1
        samIt.close();
        assertEquals(it.count, -1);
    }
}
