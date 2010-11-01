package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.BaseTest;
import static org.testng.Assert.assertEquals;
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
 * Class StingSAMIteratorTest
 * <p/>
 * Tests the StingSAMIteratorAdapter class.
 */
public class StingSAMIteratorAdapterUnitTest extends BaseTest {

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

        StingSAMIterator samIt = StingSAMIteratorAdapter.adapt(it);
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

        StingSAMIterator samIt = StingSAMIteratorAdapter.adapt(it);

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
        
        StingSAMIterator samIt = StingSAMIteratorAdapter.adapt(it);


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
