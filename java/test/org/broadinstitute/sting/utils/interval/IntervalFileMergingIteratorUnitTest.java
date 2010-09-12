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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.interval;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.Assert;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLoc;
import net.sf.picard.reference.ReferenceSequenceFileFactory;

import java.io.File;
import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Jun 14, 2010
 * Time: 10:15:52 AM
 * To change this template use File | Settings | File Templates.
 */
public class IntervalFileMergingIteratorUnitTest extends BaseTest {

        private static File refFile =  new File(validationDataLocation + "Homo_sapiens_assembly17.fasta");
        private static String intervalFileNameGATK = validationDataLocation+"test.gatk.intervals";
        private static String intervalFileNameBED = validationDataLocation+"test.bed";
        private static List<GenomeLoc> results1 = null;
        private static List<GenomeLoc> results2 = null;

        @BeforeClass
        public static void init() {
            GenomeLocParser.setupRefContigOrdering(ReferenceSequenceFileFactory.getReferenceSequenceFile(refFile));

            results1 = new ArrayList<GenomeLoc>();
            results2 = new ArrayList<GenomeLoc>();

            results1.add(GenomeLocParser.createGenomeLoc("chr1",1554));
            results1.add(GenomeLocParser.createGenomeLoc("chr1",2538,2568));
            results1.add(GenomeLocParser.createGenomeLoc("chr1",18932,19000));
            results1.add(GenomeLocParser.createGenomeLoc("chr1",19001,25000));
            results1.add(GenomeLocParser.createGenomeLoc("chr5",7415,7600));

            results2.add(GenomeLocParser.createGenomeLoc("chr1",1554));
            results2.add(GenomeLocParser.createGenomeLoc("chr1",2538,2568));
            results2.add(GenomeLocParser.createGenomeLoc("chr1",18932,25000));
            results2.add(GenomeLocParser.createGenomeLoc("chr5",7415,7600));

        }

        @Test
        public void testGATKIntervalFileIterator_Overlap() {
            logger.warn("Executing testGATKIntervalFileIterator_Overlap");

            Iterator<GenomeLoc> it = new IntervalFileMergingIterator(new File(intervalFileNameGATK),IntervalMergingRule.OVERLAPPING_ONLY);
            Iterator<GenomeLoc> check_it = results1.iterator();
            while(it.hasNext()) {
                    GenomeLoc l = it.next();
                    GenomeLoc l_expected = check_it.next();
                    //System.out.println("int: "+l+" expected: "+l_expected) ;
                    Assert.assertEquals("Unexpected location returned by the iterator: "+l,l,l_expected);
            }
       }

      @Test
      public void testGATKIntervalFileIterator_OverlapWithException() {
            logger.warn("Executing testGATKIntervalFileIterator_OverlapWithException");

            Iterator<GenomeLoc> it = new IntervalFileMergingIterator(new File(intervalFileNameGATK),IntervalMergingRule.OVERLAPPING_ONLY);
            Iterator<GenomeLoc> check_it = results1.iterator();
            try {
                while(it.hasNext()) {
                    GenomeLoc l = it.next();
                    GenomeLoc l_expected = check_it.next();
//                    System.out.println("int: "+l+" expected: "+l_expected) ;
                }
            } catch ( ReviewedStingException e) {
                    Assert.assertEquals( e.getMessage(), "Interval chr5:7414 in the interval file is out of order.");
            }
        }

     @Test
     public void testGATKIntervalFileIterator_All() {
        logger.warn("Executing testGATKIntervalFileIterator_All");

        Iterator<GenomeLoc> it = new IntervalFileMergingIterator(new File(intervalFileNameGATK),IntervalMergingRule.ALL);
        Iterator<GenomeLoc> check_it = results2.iterator();
        while(it.hasNext()) {
                GenomeLoc l = it.next();
                GenomeLoc l_expected = check_it.next();
//                System.out.println("int: "+l+" expected: "+l_expected) ;
                Assert.assertEquals("Unexpected location returned by the iterator: "+l,l,l_expected);
        }
   }

    @Test
    public void testBEDIntervalFileIterator_Overlap() {
        logger.warn("Executing testBEDIntervalFileIterator_Overlap");

        Iterator<GenomeLoc> it = new IntervalFileMergingIterator(new File(intervalFileNameBED),IntervalMergingRule.OVERLAPPING_ONLY);
        Iterator<GenomeLoc> check_it = results1.iterator();
        while(it.hasNext()) {
                GenomeLoc l = it.next();
                GenomeLoc l_expected = check_it.next();
//                System.out.println("int: "+l+" expected: "+l_expected) ;
                Assert.assertEquals("Unexpected location returned by the iterator: "+l,l,l_expected);
        }
   }

}
