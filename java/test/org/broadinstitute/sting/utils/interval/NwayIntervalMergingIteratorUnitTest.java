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

package org.broadinstitute.sting.utils.interval;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLoc;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Iterator;
import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Oct 28, 2010
 * Time: 2:46:03 PM
 * To change this template use File | Settings | File Templates.
 */
public class NwayIntervalMergingIteratorUnitTest extends BaseTest {

    private static File refFile =  new File(validationDataLocation + "Homo_sapiens_assembly17.fasta");

    private static List<GenomeLoc> stream1 = null;
    private static List<GenomeLoc> stream2 = null;
    private static List<GenomeLoc> expected = null;

    @BeforeClass
    public static void init() {
        GenomeLocParser.setupRefContigOrdering(ReferenceSequenceFileFactory.getReferenceSequenceFile(refFile));

        stream1 = new ArrayList<GenomeLoc>();
        stream2 = new ArrayList<GenomeLoc>();
        expected = new ArrayList<GenomeLoc>();

        stream1.add(GenomeLocParser.createGenomeLoc("chr1",1554,1560)); // 1
        stream1.add(GenomeLocParser.createGenomeLoc("chr1",2538,2568)); // 3
        stream1.add(GenomeLocParser.createGenomeLoc("chr1",2600,2610));    // 4
        stream1.add(GenomeLocParser.createGenomeLoc("chr1",2609,2625));    // 4
        stream1.add(GenomeLocParser.createGenomeLoc("chr1",18932,19000));  // 6
        stream1.add(GenomeLocParser.createGenomeLoc("chr1",19001,25000));  //6

        stream2.add(GenomeLocParser.createGenomeLoc("chr1",1565,1570));    //2
        stream2.add(GenomeLocParser.createGenomeLoc("chr1",2598,2604));    // 4
        stream2.add(GenomeLocParser.createGenomeLoc("chr1",7415,7600));    // 5
        stream2.add(GenomeLocParser.createGenomeLoc("chr1",18932,25000));  // 6
        stream2.add(GenomeLocParser.createGenomeLoc("chr1",30000,35000));  // 7

        expected.add(GenomeLocParser.createGenomeLoc("chr1",1554,1560)); // 1
        expected.add(GenomeLocParser.createGenomeLoc("chr1",1565,1570));    //2
        expected.add(GenomeLocParser.createGenomeLoc("chr1",2538,2568)); // 3
        expected.add(GenomeLocParser.createGenomeLoc("chr1",2598,2625));    // 4
        expected.add(GenomeLocParser.createGenomeLoc("chr1",7415,7600));    // 5
        expected.add(GenomeLocParser.createGenomeLoc("chr1",18932,25000));  // 6
        expected.add(GenomeLocParser.createGenomeLoc("chr1",30000,35000));  // 7


    }

    @Test
    public void testNwayIntervalMergingIterator() {
        logger.warn("testNwayIntervalMergingIterator");

        Iterator<GenomeLoc> it1 = stream1.iterator();
        Iterator<GenomeLoc> it2 = stream2.iterator();

        Iterator<GenomeLoc> e_it = expected.iterator();



        NwayIntervalMergingIterator it = new NwayIntervalMergingIterator(IntervalMergingRule.OVERLAPPING_ONLY);
        it.add(it1);
        it.add(it2);
        
        while(it.hasNext()) {
                GenomeLoc l = it.next();
                GenomeLoc l_expected = e_it.next();
                //System.out.println("int: "+l+" expected: "+l_expected) ;
                Assert.assertEquals(l,l_expected,"Unexpected location returned by the iterator: "+l);
        }
   }


}
