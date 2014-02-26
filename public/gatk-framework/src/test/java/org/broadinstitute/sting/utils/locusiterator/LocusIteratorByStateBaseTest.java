/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils.locusiterator;

import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.downsampling.DownsamplingMethod;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.gatk.iterators.ReadTransformer;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;

import java.util.*;

/**
 * testing of the new (non-legacy) version of LocusIteratorByState
 */
public class LocusIteratorByStateBaseTest extends BaseTest {
    protected static SAMFileHeader header;
    protected GenomeLocParser genomeLocParser;

    @BeforeClass
    public void beforeClass() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    }

    protected LocusIteratorByState makeLTBS(List<GATKSAMRecord> reads,
                                            ReadProperties readAttributes) {
        return new LocusIteratorByState(new FakeCloseableIterator<GATKSAMRecord>(reads.iterator()),
                readAttributes,
                genomeLocParser,
                LocusIteratorByState.sampleListForSAMWithoutReadGroups());
    }

    public static ReadProperties createTestReadProperties() {
        return createTestReadProperties(null, false);
    }

    public static ReadProperties createTestReadProperties( DownsamplingMethod downsamplingMethod, final boolean keepReads ) {
        return new ReadProperties(
                Collections.<SAMReaderID>emptyList(),
                new SAMFileHeader(),
                SAMFileHeader.SortOrder.coordinate,
                false,
                SAMFileReader.ValidationStringency.STRICT,
                downsamplingMethod,
                new ValidationExclusion(),
                Collections.<ReadFilter>emptyList(),
                Collections.<ReadTransformer>emptyList(),
                true,
                (byte) -1,
                keepReads);
    }

    public static class FakeCloseableIterator<T> implements CloseableIterator<T> {
        Iterator<T> iterator;

        public FakeCloseableIterator(Iterator<T> it) {
            iterator = it;
        }

        @Override
        public void close() {}

        @Override
        public boolean hasNext() {
            return iterator.hasNext();
        }

        @Override
        public T next() {
            return iterator.next();
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("Don't remove!");
        }
    }

    protected static class LIBSTest {
        public static final int locus = 44367788;
        final String cigarString;
        final int readLength;
        final private List<CigarElement> elements;

        public LIBSTest(final String cigarString) {
            final Cigar cigar = TextCigarCodec.getSingleton().decode(cigarString);
            this.cigarString = cigarString;
            this.elements = cigar.getCigarElements();
            this.readLength = cigar.getReadLength();
        }

        @Override
        public String toString() {
            return "LIBSTest{" +
                    "cigar='" + cigarString + '\'' +
                    ", readLength=" + readLength +
                    '}';
        }

        public List<CigarElement> getElements() {
            return elements;
        }

        public GATKSAMRecord makeRead() {
            GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "read", 0, locus, readLength);
            read.setReadBases(Utils.dupBytes((byte) 'A', readLength));
            final byte[] quals = new byte[readLength];
            for ( int i = 0; i < readLength; i++ )
                quals[i] = (byte)(i % QualityUtils.MAX_SAM_QUAL_SCORE);
            read.setBaseQualities(quals);
            read.setCigarString(cigarString);
            return read;
        }
    }

    private boolean isIndel(final CigarElement ce) {
        return ce.getOperator() == CigarOperator.D || ce.getOperator() == CigarOperator.I;
    }

    private boolean startsWithDeletion(final List<CigarElement> elements) {
        for ( final CigarElement element : elements ) {
            switch ( element.getOperator() ) {
                case M:
                case I:
                case EQ:
                case X:
                    return false;
                case D:
                    return true;
                default:
                    // keep looking
            }
        }

        return false;
    }

    private LIBSTest makePermutationTest(final List<CigarElement> elements) {
        CigarElement last = null;
        boolean hasMatch = false;

        // starts with D => bad
        if ( startsWithDeletion(elements) )
            return null;

        // ends with D => bad
        if ( elements.get(elements.size()-1).getOperator() == CigarOperator.D )
            return null;

        // make sure it's valid
        String cigar = "";
        int len = 0;
        for ( final CigarElement ce : elements ) {
            if ( ce.getOperator() == CigarOperator.N )
                return null; // TODO -- don't support N

            // abort on a bad cigar
            if ( last != null ) {
                if ( ce.getOperator() == last.getOperator() )
                    return null;
                if ( isIndel(ce) && isIndel(last) )
                    return null;
            }

            cigar += ce.getLength() + ce.getOperator().toString();
            len += ce.getLength();
            last = ce;
            hasMatch = hasMatch || ce.getOperator() == CigarOperator.M;
        }

        if ( ! hasMatch && elements.size() == 1 &&
                ! (last.getOperator() == CigarOperator.I || last.getOperator() == CigarOperator.S))
            return null;

        return new LIBSTest(cigar);
    }

    @DataProvider(name = "LIBSTest")
    public Object[][] createLIBSTests(final List<Integer> cigarLengths, final List<Integer> combinations) {
        final List<Object[]> tests = new LinkedList<Object[]>();

        final List<CigarOperator> allOps = Arrays.asList(CigarOperator.values());

        final List<CigarElement> singleCigars = new LinkedList<CigarElement>();
        for ( final int len : cigarLengths )
            for ( final CigarOperator op : allOps )
                singleCigars.add(new CigarElement(len, op));

        for ( final int complexity : combinations ) {
            for ( final List<CigarElement> elements : Utils.makePermutations(singleCigars, complexity, true) ) {
                final LIBSTest test = makePermutationTest(elements);
                if ( test != null ) tests.add(new Object[]{test});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    /**
     * Work around inadequate tests that aren't worth fixing.
     *
     * Look at the CIGAR 2M2P2D2P2M.  Both M states border a deletion, separated by P (padding elements).  So
     * the right answer for deletions here is true for isBeforeDeletion() and isAfterDeletion() for the first
     * and second M.  But the LIBS_position doesn't say so.
     *
     * @param elements
     * @return
     */
    protected static boolean hasNeighboringPaddedOps(final List<CigarElement> elements, final int elementI) {
        return (elementI - 1 >= 0 && isPadding(elements.get(elementI-1))) ||
                (elementI + 1 < elements.size() && isPadding(elements.get(elementI+1)));
    }

    private static boolean isPadding(final CigarElement elt) {
        return elt.getOperator() == CigarOperator.P || elt.getOperator() == CigarOperator.H || elt.getOperator() == CigarOperator.S;
    }
}
