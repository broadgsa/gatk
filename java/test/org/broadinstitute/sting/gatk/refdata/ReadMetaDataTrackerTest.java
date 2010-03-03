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
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashSet;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.TreeMap;


/**
 * @author aaron
 *         <p/>
 *         Class ReadMetaDataTrackerTest
 *         <p/>
 *         test out the ReadMetaDataTracker
 */
public class ReadMetaDataTrackerTest extends BaseTest {
    private static int startingChr = 1;
    private static int endingChr = 2;
    private static int readCount = 100;
    private static int DEFAULT_READ_LENGTH = ArtificialSAMUtils.DEFAULT_READ_LENGTH;
    private static SAMFileHeader header;

    @BeforeClass
    public static void beforeClass() {
        header = ArtificialSAMUtils.createArtificialSamHeader((endingChr - startingChr) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);
        GenomeLocParser.setupRefContigOrdering(header.getSequenceDictionary());
    }

    @Before
    public void beforeEach() {
    }

    @Test
    public void rodAtEachReadBase() {
       ReadMetaDataTracker tracker = getRMDT(1);

        // count the positions
        int count = 0;
        for (int x : tracker.getPositionMapping().keySet()) {
            count++;
            Assert.assertEquals(1, tracker.getPositionMapping().get(x).size());
        }
        Assert.assertEquals(10, count);
    }

    @Test
    public void sparceRODsForRead() {
        ReadMetaDataTracker tracker = getRMDT(7);

        // count the positions
        int count = 0;
        for (int x : tracker.getPositionMapping().keySet()) {
            count++;
            Assert.assertEquals(1, tracker.getPositionMapping().get(x).size());
        }
        Assert.assertEquals(2, count);
    }

    @Test
    public void rodByGenomeLoc() {
        ReadMetaDataTracker tracker = getRMDT(1);

        // count the positions
        int count = 0;
        for (Long x : tracker.getGenomeLocMapping().keySet()) {
            count++;
            Assert.assertEquals(1, tracker.getGenomeLocMapping().get(x).size());
        }
        Assert.assertEquals(10, count);
    }

    private ReadMetaDataTracker getRMDT(int incr) {
        SAMRecord record = ArtificialSAMUtils.createArtificialRead(header, "name", 0, 1, 10);
        byte[] c = new byte[10];
        for (int x = 0; x < 10; x++)
            c[x] = 'A';
        record.setReadBases(c);
        TreeMap<Long, Set<ReferenceOrderedDatum>> data = new TreeMap<Long, Set<ReferenceOrderedDatum>>();
        for (int x = 0; x < record.getAlignmentEnd(); x+=incr) {
            GenomeLoc loc = GenomeLocParser.createGenomeLoc(record.getReferenceIndex(), record.getAlignmentStart() + x, record.getAlignmentStart() + x);
            Set<ReferenceOrderedDatum> set = new HashSet<ReferenceOrderedDatum>();
            set.add(new FakeRODatum(loc));
            data.put((long)record.getAlignmentStart() + x,set);
        }
        ReadMetaDataTracker tracker = new ReadMetaDataTracker(record, data);
        return tracker;
    }


    /** for testing only */
    static public class FakeRODatum implements ReferenceOrderedDatum {

        final GenomeLoc location;

        public FakeRODatum(GenomeLoc location) {
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
    }
}