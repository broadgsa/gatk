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

package org.broadinstitute.gatk.engine.filters;


import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.utils.ValidationExclusion;
import org.broadinstitute.gatk.engine.datasources.reads.SAMDataSource;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.exceptions.UserException.UnsupportedCigarOperatorException;

import java.lang.annotation.*;
import java.lang.reflect.Method;
import java.util.*;


/**
 * Tests for the MalformedReadFilter
 *
 * @author Eric Banks
 * @since 3/14/13
 */
public class MalformedReadFilterUnitTest extends ReadFilterTest {

    //////////////////////////////////////
    // Test the checkSeqStored() method //
    //////////////////////////////////////

    @Test(enabled = true)
    public void testCheckSeqStored () {

        final GATKSAMRecord goodRead = ArtificialSAMUtils.createArtificialRead(new byte[]{(byte)'A'}, new byte[]{(byte)'A'}, "1M");
        final GATKSAMRecord badRead = ArtificialSAMUtils.createArtificialRead(new byte[]{}, new byte[]{}, "1M");
        badRead.setReadString("*");

        Assert.assertTrue(MalformedReadFilter.checkSeqStored(goodRead, true));
        Assert.assertFalse(MalformedReadFilter.checkSeqStored(badRead, true));

        try {
            MalformedReadFilter.checkSeqStored(badRead, false);
            Assert.assertTrue(false, "We should have exceptioned out in the previous line");
        } catch (UserException e) { }
    }

    @Test(enabled = true, dataProvider= "UnsupportedCigarOperatorDataProvider")
    @CigarOperatorTest(CigarOperatorTest.Outcome.FILTER)
    public void testCigarNOperatorFilterTruePositive(String cigarString) {

       final MalformedReadFilter filter = buildMalformedReadFilter(true);
       final SAMRecord nContainingCigarRead = buildSAMRecord(cigarString);
       Assert.assertTrue(filter.filterOut(nContainingCigarRead),
                  " Did not filtered out a N containing CIGAR read");
    }

    @Test(enabled = true, dataProvider= "UnsupportedCigarOperatorDataProvider")
    @CigarOperatorTest(CigarOperatorTest.Outcome.ACCEPT)
    public void testCigarNOperatorFilterTrueNegative(String cigarString) {

        final MalformedReadFilter filter = buildMalformedReadFilter(true);
        final SAMRecord nonNContainingCigarRead = buildSAMRecord(cigarString);
        Assert.assertFalse(filter.filterOut(nonNContainingCigarRead),
                    " Filtered out a non-N containing CIGAR read");
    }

    @Test(enabled = true,
            expectedExceptions = UnsupportedCigarOperatorException.class,
            dataProvider= "UnsupportedCigarOperatorDataProvider")
    @CigarOperatorTest(CigarOperatorTest.Outcome.EXCEPTION)
    public void testCigarNOperatorFilterException(final String cigarString) {

        final MalformedReadFilter filter = buildMalformedReadFilter(false);
        final SAMRecord nContainingCigarRead = buildSAMRecord(cigarString);

        filter.filterOut(nContainingCigarRead);
    }

    @Test(enabled = true, dataProvider="UnsupportedCigarOperatorDataProvider")
    @CigarOperatorTest(CigarOperatorTest.Outcome.ACCEPT)
    public void testCigarNOperatorFilterControl(final String cigarString) {

        final MalformedReadFilter filter = buildMalformedReadFilter(false);
        final SAMRecord nonNContainingCigarRead = buildSAMRecord(cigarString);

        Assert.assertFalse(filter.filterOut(nonNContainingCigarRead));
    }

    protected SAMRecord buildSAMRecord(final String cigarString) {
        final Cigar nContainingCigar = TextCigarCodec.decode(cigarString);
        return  this.createRead(nContainingCigar, 1, 0, 10);
    }

    protected MalformedReadFilter buildMalformedReadFilter(final boolean filterRNO) {
        return buildMalformedReadFiter(filterRNO,new ValidationExclusion.TYPE[] {});
    }

    protected MalformedReadFilter buildMalformedReadFiter(boolean filterRNO, final ValidationExclusion.TYPE... excl) {
        final ValidationExclusion ve = new ValidationExclusion(Arrays.asList(excl));

        final MalformedReadFilter filter = new MalformedReadFilter();

        final SAMFileHeader h = getHeader();
        final SAMDataSource ds =  getDataSource();

        final GenomeAnalysisEngine gae = new GenomeAnalysisEngine() {
            @Override
            public SAMFileHeader getSAMFileHeader() {
                return h;
            }

            @Override
            public SAMDataSource getReadsDataSource() {
                return ds;
            }
        };
        filter.initialize(gae);
        filter.filterReadsWithNCigar = filterRNO;
        return filter;
    }

    @Retention(RetentionPolicy.RUNTIME)
    @Target(ElementType.METHOD)
    @Inherited
    protected @interface CigarOperatorTest {

        enum Outcome {
            ANY,ACCEPT,FILTER,EXCEPTION,IGNORE;

            public boolean appliesTo (String cigar) {
                boolean hasN = cigar.indexOf('N') != -1;
                switch (this) {
                    case ANY: return true;
                    case ACCEPT: return !hasN;
                    case IGNORE: return hasN;
                    case FILTER:
                    case EXCEPTION:
                    default:
                        return hasN;

                }
            }
        }

        Outcome value() default Outcome.ANY;
    }

    /**
     * Cigar test data for unsupported operator test.
     * Each element of this array corresponds to a test case. In turn the first element of the test case array is the
     * Cigar string for that test case and the second indicates whether it should be filtered due to the presence of a
     * unsupported operator
     */
    private static final String[] TEST_CIGARS =  {
       "101M10D20I10M",
       "6M14N5M",
       "1N",
       "101M",
       "110N",
       "2N4M",
       "4M2N",
       "3M1I1M",
       "1M2I2M",
       "1M10N1I1M",
       "1M1I1D",
       "11N12M1I34M12N"
    };

    @DataProvider(name= "UnsupportedCigarOperatorDataProvider")
    public Iterator<Object[]> unsupportedOperatorDataProvider(final Method testMethod) {
        final CigarOperatorTest a = resolveCigarOperatorTestAnnotation(testMethod);
        final List<Object[]> result = new LinkedList<Object[]>();
        for (final String cigarString : TEST_CIGARS) {
            if (a == null || a.value().appliesTo(cigarString)) {
                result.add(new Object[] { cigarString });
            }
        }
        return result.iterator();
    }

    /**
     * Gets the most specific {@link CigarOperatorTest} annotation for the
     * signature of the test method provided.
     * <p/>
     * This in-house implementation is required due to the fact that method
     * annotations do not have inheritance.
     *
     * @param m targeted test method.
     * @return <code>null</code> if there is no {@link CigarOperatorTest}
     * annotation in this or overridden methods.
     */
    private CigarOperatorTest resolveCigarOperatorTestAnnotation(final Method m) {
       CigarOperatorTest res = m.getAnnotation(CigarOperatorTest.class);
       if (res != null) {
           return res;
       }
       Class<?> c = this.getClass();
       Class<?> p = c.getSuperclass();
       while (p != null && p != Object.class) {
           try {
             final Method met = p.getDeclaredMethod(m.getName(),
                     m.getParameterTypes());
             res = met.getAnnotation(CigarOperatorTest.class);
             if (res != null) {
                 break;
             }
           } catch (NoSuchMethodException e) {
             // Its ok; nothing to do here, just keep looking.
           }
           c = p;
           p = c.getSuperclass();
       }
       return res;
    }

}
