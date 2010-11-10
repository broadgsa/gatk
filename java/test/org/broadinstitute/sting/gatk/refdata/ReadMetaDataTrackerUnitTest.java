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

package org.broadinstitute.sting.gatk.refdata;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import org.testng.Assert;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.datasources.providers.RODMetaDataContainer;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
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
 *         Class ReadMetaDataTrackerUnitTest
 *         <p/>
 *         test out the ReadMetaDataTracker
 */
public class ReadMetaDataTrackerUnitTest extends BaseTest {
    private static int startingChr = 1;
    private static int endingChr = 2;
    private static int readCount = 100;
    private static int DEFAULT_READ_LENGTH = ArtificialSAMUtils.DEFAULT_READ_LENGTH;
    private static SAMFileHeader header;
    private Set<String> nameSet;

    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void beforeClass() {
        header = ArtificialSAMUtils.createArtificialSamHeader((endingChr - startingChr) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    }

    @BeforeMethod
    public void beforeEach() {
        nameSet = new TreeSet<String>();
        nameSet.add("default");
    }

    @Test
    public void twoRodsAtEachReadBase() {
        nameSet.add("default2");
        ReadMetaDataTracker tracker = getRMDT(1, nameSet, true);

        // count the positions
        int count = 0;
        for (Integer x : tracker.getReadOffsetMapping().keySet()) {
            count++;
            Assert.assertEquals(tracker.getReadOffsetMapping().get(x).size(), 2);
        }
        Assert.assertEquals(count, 10);
    }

    @Test
    public void rodAtEachReadBase() {

        ReadMetaDataTracker tracker = getRMDT(1, nameSet, true);

        // count the positions
        int count = 0;
        for (Integer x : tracker.getReadOffsetMapping().keySet()) {
            count++;
            Assert.assertEquals(tracker.getReadOffsetMapping().get(x).size(), 1);
        }
        Assert.assertEquals(count, 10);
    }

    @Test
    public void filterByName() {
        nameSet.add("default2");
        ReadMetaDataTracker tracker = getRMDT(1, nameSet, true);

        // count the positions
        int count = 0;
        Map<Integer, Collection<GATKFeature>> map = tracker.getReadOffsetMapping("default");
        for (Integer x : map.keySet()) {
            count++;
            Assert.assertEquals(map.get(x).size(), 1);
        }
        Assert.assertEquals(count, 10);
    }

    @Test
    public void filterByDupType() {
        nameSet.add("default2");
        ReadMetaDataTracker tracker = getRMDT(1, nameSet, false);  // create both RODs of the same type
        // count the positions
        int count = 0;
        Map<Integer, Collection<GATKFeature>> map = tracker.getReadOffsetMapping(FakeRODatum.class);
        for (Integer x : map.keySet()) {
            count++;
            Assert.assertEquals(map.get(x).size(), 2);
        }
        Assert.assertEquals(count, 10);
    }

    // @Test this test can be uncommented to determine the speed impacts of any changes to the RODs for reads system

    public void filterByMassiveDupType() {

        for (int y = 0; y < 20; y++) {
            nameSet.add("default" + String.valueOf(y));
            long firstTime = System.currentTimeMillis();
            for (int lp = 0; lp < 1000; lp++) {
                ReadMetaDataTracker tracker = getRMDT(1, nameSet, false);  // create both RODs of the same type
                // count the positions
                int count = 0;
                Map<Integer, Collection<GATKFeature>> map = tracker.getReadOffsetMapping(FakeRODatum.class);
                for (Integer x : map.keySet()) {
                    count++;
                    Assert.assertEquals(map.get(x).size(), y + 2);
                }
                Assert.assertEquals(count, 10);
            }
            System.err.println(y + " = " + (System.currentTimeMillis() - firstTime));
        }
    }


    @Test
    public void filterByType() {
        nameSet.add("default2");
        ReadMetaDataTracker tracker = getRMDT(1, nameSet, true);

        // count the positions
        int count = 0;
        Map<Integer, Collection<GATKFeature>> map = tracker.getReadOffsetMapping(Fake2RODatum.class);
        for (int x : map.keySet()) {
            count++;
            Assert.assertEquals(map.get(x).size(), 1);
        }
        Assert.assertEquals(count, 10);
    }

    @Test
    public void sparceRODsForRead() {
        ReadMetaDataTracker tracker = getRMDT(7, nameSet, true);

        // count the positions
        int count = 0;
        for (Integer x : tracker.getReadOffsetMapping().keySet()) {
            count++;
            Assert.assertEquals(tracker.getReadOffsetMapping().get(x).size(), 1);
        }
        Assert.assertEquals(count, 2);
    }

    @Test
    public void rodByGenomeLoc() {
        ReadMetaDataTracker tracker = getRMDT(1, nameSet, true);

        // count the positions
        int count = 0;
        for (Integer x : tracker.getContigOffsetMapping().keySet()) {
            count++;
            Assert.assertEquals(tracker.getContigOffsetMapping().get(x).size(), 1);
        }
        Assert.assertEquals(count, 10);
    }


    /**
     * create a ReadMetaDataTracker given:
     *
     * @param incr  the spacing between site locations
     * @param names the names of the reference ordered data to create: one will be created at every location for each name
     *
     * @return a ReadMetaDataTracker
     */
    private ReadMetaDataTracker getRMDT(int incr, Set<String> names, boolean alternateTypes) {
        SAMRecord record = ArtificialSAMUtils.createArtificialRead(header, "name", 0, 1, 10);
        TreeMap<Integer, RODMetaDataContainer> data = new TreeMap<Integer, RODMetaDataContainer>();
        for (int x = 0; x < record.getAlignmentEnd(); x += incr) {
            GenomeLoc loc = genomeLocParser.createGenomeLoc(record.getReferenceName(), record.getAlignmentStart() + x, record.getAlignmentStart() + x);
            RODMetaDataContainer set = new RODMetaDataContainer();

            int cnt = 0;
            for (String name : names) {
                if (alternateTypes)
                    set.addEntry((cnt % 2 == 0) ? new FakeRODatum(loc, name) : new Fake2RODatum(loc, name));
                else
                    set.addEntry(new FakeRODatum(loc, name));
                cnt++;
            }
            data.put(record.getAlignmentStart() + x, set);
        }
        ReadMetaDataTracker tracker = new ReadMetaDataTracker(genomeLocParser, record, data);
        return tracker;
    }


    /** for testing, we want a fake rod with a different classname, for the get-by-class-name functions */
    static public class Fake2RODatum extends FakeRODatum {

        public Fake2RODatum(GenomeLoc location, String name) {
            super(location, name);
        }
    }


    /** for testing only */
    static public class FakeRODatum extends GATKFeature {

        final GenomeLoc location;
        final String name;

        public FakeRODatum(GenomeLoc location, String name) {
            super(name);
            this.location = location;
            this.name = name;
        }

        @Override
        public String getName() {
            return name;
        }

        @Override
        public GenomeLoc getLocation() {
            return this.location;
        }

        @Override
        public Object getUnderlyingObject() {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public String getChr() {
            return location.getContig();
        }

        @Override
        public int getStart() {
            return (int)this.location.getStart();
        }

        @Override
        public int getEnd() {
            return (int)this.location.getStop();
        }
    }
}
