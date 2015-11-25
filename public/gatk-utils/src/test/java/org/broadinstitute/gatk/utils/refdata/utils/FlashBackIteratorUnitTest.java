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

package org.broadinstitute.gatk.utils.refdata.utils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.testng.Assert;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.refdata.ReferenceOrderedDatum;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;

import org.testng.annotations.BeforeMethod;

import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.List;


/**
 * @author aaron
 *         <p/>
 *         Class FlashBackIteratorUnitTest
 *         <p/>
 *         just like a greatful dead show...this will be prone to flashbacks
 */
public class FlashBackIteratorUnitTest extends BaseTest {
    private SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(NUMBER_OF_CHROMOSOMES, STARTING_CHROMOSOME, CHROMOSOME_SIZE);
    private static final int NUMBER_OF_CHROMOSOMES = 5;
    private static final int STARTING_CHROMOSOME = 1;
    private static final int CHROMOSOME_SIZE = 1000;

    private String firstContig;
    private GenomeLocParser genomeLocParser;

    @BeforeMethod
    public void setup() {
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
        firstContig = header.getSequenceDictionary().getSequence(0).getSequenceName();
    }

    @Test
    public void testBasicIteration() {
        GenomeLoc loc = genomeLocParser.createGenomeLoc(firstContig, 0, 0);
        FlashBackIterator iter = new FlashBackIterator(new FakeSeekableRODIterator(genomeLocParser,loc));
        GenomeLoc lastLocation = null;
        for (int x = 0; x < 10; x++) {
            iter.next();
            GenomeLoc cur = iter.position();
            if (lastLocation != null) {
                Assert.assertTrue(lastLocation.isBefore(cur));
            }
            lastLocation = cur;
        }
    }

    @Test
    public void testBasicIterationThenFlashBack() {
        GenomeLoc loc = genomeLocParser.createGenomeLoc(firstContig, 0, 0);
        FlashBackIterator iter = new FlashBackIterator(new FakeSeekableRODIterator(genomeLocParser,loc));
        GenomeLoc lastLocation = null;
        for (int x = 0; x < 10; x++) {
            iter.next();
            GenomeLoc cur = iter.position();
            if (lastLocation != null) {
                Assert.assertTrue(lastLocation.isBefore(cur));
            }
            lastLocation = cur;
        }
        iter.flashBackTo(genomeLocParser.createGenomeLoc(firstContig, 2));
    }

    @Test
    public void testBasicIterationThenFlashBackThenIterate() {
        GenomeLoc loc = genomeLocParser.createGenomeLoc(firstContig, 0, 0);
        FlashBackIterator iter = new FlashBackIterator(new FakeSeekableRODIterator(genomeLocParser,loc));
        GenomeLoc lastLocation = null;
        for (int x = 0; x < 10; x++) {
            iter.next();
            GenomeLoc cur = iter.position();
            if (lastLocation != null) {
                Assert.assertTrue(lastLocation.isBefore(cur));
            }
            lastLocation = cur;
        }
        iter.flashBackTo(genomeLocParser.createGenomeLoc(firstContig, 1));
        int count = 0;
        while (iter.hasNext()) {
            count++;
            iter.next();
        }
        Assert.assertEquals(count, 10);
    }


    @Test
    public void testFlashBackTruth() {
        GenomeLoc loc = genomeLocParser.createGenomeLoc(firstContig, 0, 0);
        LocationAwareSeekableRODIterator backIter = new FakeSeekableRODIterator(genomeLocParser,loc);
        // remove the first three records
        backIter.next();
        backIter.next();
        backIter.next();
        FlashBackIterator iter = new FlashBackIterator(backIter);
        GenomeLoc lastLocation = null;
        for (int x = 0; x < 10; x++) {
            iter.next();
            GenomeLoc cur = iter.position();
            if (lastLocation != null) {
                Assert.assertTrue(lastLocation.isBefore(cur));
            }
            lastLocation = cur;
        }
        Assert.assertTrue(iter.canFlashBackTo(genomeLocParser.createGenomeLoc(firstContig, 5)));
        Assert.assertTrue(iter.canFlashBackTo(genomeLocParser.createGenomeLoc(firstContig, 15)));
        Assert.assertTrue(!iter.canFlashBackTo(genomeLocParser.createGenomeLoc(firstContig, 2)));
        Assert.assertTrue(!iter.canFlashBackTo(genomeLocParser.createGenomeLoc(firstContig, 1)));
    }

    @Test
    public void testBasicIterationThenFlashBackHalfWayThenIterate() {
        GenomeLoc loc = genomeLocParser.createGenomeLoc(firstContig, 0, 0);
        FlashBackIterator iter = new FlashBackIterator(new FakeSeekableRODIterator(genomeLocParser,loc));
        GenomeLoc lastLocation = null;
        for (int x = 0; x < 10; x++) {
            iter.next();
            GenomeLoc cur = iter.position();
            if (lastLocation != null) {
                Assert.assertTrue(lastLocation.isBefore(cur));
            }
            lastLocation = cur;
        }
        iter.flashBackTo(genomeLocParser.createGenomeLoc(firstContig, 5));
        int count = 0;
        while (iter.hasNext()) {
            count++;
            iter.next();
        }
        Assert.assertEquals(count, 6); // chr1:5, 6, 7, 8, 9, and 10
    }
}


class FakeSeekableRODIterator implements LocationAwareSeekableRODIterator {
    private GenomeLocParser genomeLocParser;

    // current location
    private GenomeLoc location;
    private FakeRODatum curROD;
    private int recordCount = 10;

    public FakeSeekableRODIterator(GenomeLocParser genomeLocParser,GenomeLoc startingLoc) {
        this.genomeLocParser = genomeLocParser;
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
        this.location = interval;
        return next();
    }

    @Override
    public boolean hasNext() {
        return (recordCount > 0);
    }

    @Override
    public RODRecordList next() {
        RODRecordList list = new FakeRODRecordList();
        curROD = new FakeRODatum("STUPIDNAME", location);
        location = genomeLocParser.createGenomeLoc(location.getContig(), location.getStart() + 1, location.getStop() + 1);
        list.add(curROD);
        recordCount--;
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


/** for testing only */
class FakeRODatum extends GATKFeature implements ReferenceOrderedDatum {

    final GenomeLoc location;

    public FakeRODatum(String name, GenomeLoc location) {
        super(name);
        this.location = location;
    }

    @Override
    public String getName() {
        return "false";
    }

    @Override
    public boolean parseLine(Object header, String[] parts) throws IOException {
        return false;
    }

    @Override
    public String toSimpleString() {
        return "";
    }

    @Override
    public String repl() {
        return "";
    }

    /**
     * Used by the ROD system to determine how to split input lines
     *
     * @return Regex string delimiter separating fields
     */
    @Override
    public String delimiterRegex() {
        return "";
    }

    @Override
    public GenomeLoc getLocation() {
        return location;
    }

    @Override
    public Object getUnderlyingObject() {
        return this;
    }

    @Override
    public int compareTo(ReferenceOrderedDatum that) {
        return location.compareTo(that.getLocation());
    }

    /**
     * Backdoor hook to read header, meta-data, etc. associated with the file.  Will be
     * called by the ROD system before streaming starts
     *
     * @param source source data file on disk from which this rod stream will be pulled
     *
     * @return a header object that will be passed to parseLine command
     */
    @Override
    public Object initialize(File source) throws FileNotFoundException {
        return null;
    }

    @Override
    public String getChr() {
        return getContig();
    }

    @Override
    public String getContig() {
        return location.getContig();
    }

    @Override
    public int getStart() {
        return (int)location.getStart();
    }

    @Override
    public int getEnd() {
        return (int)location.getStop();
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
