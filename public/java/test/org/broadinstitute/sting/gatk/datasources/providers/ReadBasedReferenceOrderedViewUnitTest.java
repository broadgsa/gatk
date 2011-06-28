/*
 * Copyright (c) 2010.  The Broad Institute
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
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.datasources.providers;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import org.testng.Assert;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTrackerUnitTest;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;

import org.testng.annotations.BeforeMethod;

import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.*;


/**
 * @author aaron
 *         <p/>
 *         Class ReadBasedReferenceOrderedViewUnitTest
 *         <p/>
 *         test out the ReadBasedReferenceOrderedView class
 */
public class ReadBasedReferenceOrderedViewUnitTest extends BaseTest {
    private GenomeLocParser genomeLocParser;

    private static int startingChr = 1;
    private static int endingChr = 2;
    private static int readCount = 100;
    private static int DEFAULT_READ_LENGTH = ArtificialSAMUtils.DEFAULT_READ_LENGTH;
    private static SAMFileHeader header;

    @BeforeClass
    public void beforeClass() {
        header = ArtificialSAMUtils.createArtificialSamHeader((endingChr - startingChr) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    }

    @BeforeMethod
    public void beforeEach() {
    }

    @Test
    public void testCreateReadMetaDataTrackerOnePerSite() {
        // make ten reads,
        List<SAMRecord> records = new ArrayList<SAMRecord>();
        for (int x = 1; x < 11; x++) {
            SAMRecord rec = ArtificialSAMUtils.createArtificialRead(header, "name", 0, x, 10);
        }
        GenomeLoc start = genomeLocParser.createGenomeLoc(header.getSequenceDictionary().getSequence(0).getSequenceName(), 0, 0);
        List<RMDDataState> list = new ArrayList<RMDDataState>();
        list.add(new RMDDataState(null, new FakePeekingRODIterator(genomeLocParser,start, "fakeName")));
        ReadBasedReferenceOrderedView view = new ReadBasedReferenceOrderedView(new WindowedData(list));

        for (SAMRecord rec : records) {
            ReadMetaDataTracker tracker = view.getReferenceOrderedDataForRead(rec);
            Map<Integer, Collection<GATKFeature>> map = tracker.getReadOffsetMapping();
            for (Integer i : map.keySet()) {
                Assert.assertEquals(map.get(i).size(), 1);
            }
            Assert.assertEquals(map.keySet().size(), 10);
        }

    }

}


class FakePeekingRODIterator implements LocationAwareSeekableRODIterator {
    private GenomeLocParser genomeLocParser;

    // current location
    private GenomeLoc location;
    private GATKFeature curROD;
    private final String name;

    public FakePeekingRODIterator(GenomeLocParser genomeLocParser, GenomeLoc startingLoc, String name) {
        this.name = name;
        this.location = genomeLocParser.createGenomeLoc(startingLoc.getContig(), startingLoc.getStart() + 1, startingLoc.getStop() + 1);
    }

    /**
     * Gets the header associated with the backing input stream.
     * @return the ROD header.
     */
    @Override
    public Object getHeader() {
        return null;
    }

    /**
     * Gets the sequence dictionary associated with the backing input stream.
     * @return sequence dictionary from the ROD header.
     */
    @Override
    public SAMSequenceDictionary getSequenceDictionary() {
        return null;
    }


    @Override
    public GenomeLoc peekNextLocation() {
        System.err.println("Peek Next -> " + location);
        return location;
    }

    @Override
    public GenomeLoc position() {
        return location;
    }

    @Override
    public RODRecordList seekForward(GenomeLoc interval) {
        while (location.isBefore(interval))
            next();
        return next(); // we always move by one, we know the next location will be right
    }

    @Override
    public boolean hasNext() {
        return true; // we always have next
    }

    @Override
    public RODRecordList next() {
        System.err.println("Next -> " + location);
        curROD = new ReadMetaDataTrackerUnitTest.FakeRODatum(location, name);
        location = genomeLocParser.createGenomeLoc(location.getContig(), location.getStart() + 1, location.getStop() + 1);
        FakeRODRecordList list = new FakeRODRecordList();
        list.add(curROD);
        return list;
    }

    @Override
    public void remove() {
        throw new IllegalStateException("GRRR");
    }

    @Override
    public void close() {
        // nothing to do
    }
}

class FakeRODRecordList extends AbstractList<GATKFeature> implements RODRecordList {
    private final List<GATKFeature> list = new ArrayList<GATKFeature>();

    public boolean add(GATKFeature data) {
        return list.add(data);
    }

    @Override
    public GATKFeature get(int i) {
        return list.get(i);
    }

    @Override
    public int size() {
        return list.size();
    }

    @Override
    public GenomeLoc getLocation() {
        return list.get(0).getLocation();
    }

    @Override
    public String getName() {
        return "test";
    }

    @Override
    public int compareTo(RODRecordList rodRecordList) {
        return this.list.get(0).getLocation().compareTo(rodRecordList.getLocation());
    }
}
