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

import org.testng.Assert;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMReadGroupRecord;

import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

public class ReadGroupBlackListFilterUnitTest extends ReadFilterTest {

    @Test(expectedExceptions=ReviewedGATKException.class)
    public void testBadFilter() {
        List<String> badFilters = Collections.singletonList("bad");
        new ReadGroupBlackListFilter(badFilters);
    }
    @Test(expectedExceptions=ReviewedGATKException.class)
    public void testBadFilterTag() {
        List<String> badFilters = Collections.singletonList("bad:filter");
        new ReadGroupBlackListFilter(badFilters);
    }

    @Test(expectedExceptions=ReviewedGATKException.class)
    public void testBadFilterFile() {
        List<String> badFilters = Collections.singletonList("/foo/bar/rgbl.txt");
        new ReadGroupBlackListFilter(badFilters);
    }

    @Test
    public void testFilterReadGroup() {
        SAMRecord filteredRecord = ArtificialSAMUtils.createArtificialRead(getHeader(), "readUno", 0, 1, 20);
        filteredRecord.setAttribute("RG", getReadGroupId(1));

        SAMRecord unfilteredRecord = ArtificialSAMUtils.createArtificialRead(getHeader(), "readDos", 0, 2, 20);
        unfilteredRecord.setAttribute("RG", getReadGroupId(2));

        List<String> filterList = new ArrayList<String>();
        filterList.add("RG:" + getReadGroupId(1));

        ReadGroupBlackListFilter filter = new ReadGroupBlackListFilter(filterList);
        Assert.assertTrue(filter.filterOut(filteredRecord));
        Assert.assertFalse(filter.filterOut(unfilteredRecord));
    }

    @Test
    public void testFilterPlatformUnit() {
        SAMRecord filteredRecord = ArtificialSAMUtils.createArtificialRead(getHeader(), "readUno", 0, 1, 20);
        filteredRecord.setAttribute("RG", getReadGroupId(1));

        SAMRecord unfilteredRecord = ArtificialSAMUtils.createArtificialRead(getHeader(), "readDos", 0, 2, 20);
        unfilteredRecord.setAttribute("RG", getReadGroupId(2));

        List<String> filterList = new ArrayList<String>();
        filterList.add("PU:" + getPlatformUnit(1));

        ReadGroupBlackListFilter filter = new ReadGroupBlackListFilter(filterList);
        Assert.assertTrue(filter.filterOut(filteredRecord));
        Assert.assertFalse(filter.filterOut(unfilteredRecord));
    }

    @Test
    public void testFilterOutByReadGroup() {
        int recordsPerGroup = 3;
        List<SAMRecord> records = new ArrayList<SAMRecord>();
        int alignmentStart = 0;
        for (int x = 1; x <= getReadGroupCount(); x++) {
            SAMReadGroupRecord groupRecord = getHeader().getReadGroup(getReadGroupId(x));
            for (int y = 1; y <= recordsPerGroup; y++) {
                SAMRecord record = ArtificialSAMUtils.createArtificialRead(getHeader(), "readUno", 0, ++alignmentStart, 20);
                record.setAttribute("RG", groupRecord.getReadGroupId());
                records.add(record);
            }
        }

        List<String> filterList = new ArrayList<String>();
        filterList.add("RG:" + getReadGroupId(1));
        filterList.add("RG:" + getReadGroupId(3));

        ReadGroupBlackListFilter filter = new ReadGroupBlackListFilter(filterList);
        int filtered = 0;
        int unfiltered = 0;
        for (SAMRecord record : records) {
            String readGroupName = record.getReadGroup().getReadGroupId();
            if (filter.filterOut(record)) {
                if (!filterList.contains("RG:" + readGroupName))
                    Assert.fail("Read group " + readGroupName + " was filtered");
                filtered++;
            } else {
                if (filterList.contains("RG:" + readGroupName))
                    Assert.fail("Read group " + readGroupName + " was not filtered");
                unfiltered++;
            }
        }

        int filteredExpected = recordsPerGroup * 2;
        int unfilteredExpected = recordsPerGroup * (getReadGroupCount() - 2);
        Assert.assertEquals(filtered, filteredExpected, "Filtered");
        Assert.assertEquals(unfiltered, unfilteredExpected, "Uniltered");
    }

    @Test
    public void testFilterOutByAttribute() {
        int recordsPerGroup = 3;
        List<SAMRecord> records = new ArrayList<SAMRecord>();
        int alignmentStart = 0;
        for (int x = 1; x <= getReadGroupCount(); x++) {
            SAMReadGroupRecord groupRecord = getHeader().getReadGroup(getReadGroupId(x));
            for (int y = 1; y <= recordsPerGroup; y++) {
                SAMRecord record = ArtificialSAMUtils.createArtificialRead(getHeader(), "readUno", 0, ++alignmentStart, 20);
                record.setAttribute("RG", groupRecord.getReadGroupId());
                records.add(record);
            }
        }

        List<String> filterList = new ArrayList<String>();
        filterList.add("PU:" + getPlatformUnit(1));

        ReadGroupBlackListFilter filter = new ReadGroupBlackListFilter(filterList);
        int filtered = 0;
        int unfiltered = 0;
        for (SAMRecord record : records) {
            String platformUnit = (String) record.getReadGroup().getAttribute("PU");
            if (filter.filterOut(record)) {
                if (!filterList.contains("PU:" + platformUnit))
                    Assert.fail("Platform unit " + platformUnit + " was filtered");
                filtered++;
            } else {
                if (filterList.contains("PU:" + platformUnit))
                    Assert.fail("Platform unit " + platformUnit + " was not filtered");
                unfiltered++;
            }
        }

        int filteredExpected = 6;
        int unfilteredExpected = 9;
        Assert.assertEquals(filtered, filteredExpected, "Filtered");
        Assert.assertEquals(unfiltered, unfilteredExpected, "Uniltered");
    }

    @Test
    public void testFilterOutByFile() {
        int recordsPerGroup = 3;
        List<SAMRecord> records = new ArrayList<SAMRecord>();
        int alignmentStart = 0;
        for (int x = 1; x <= getReadGroupCount(); x++) {
            SAMReadGroupRecord groupRecord = getHeader().getReadGroup(getReadGroupId(x));
            for (int y = 1; y <= recordsPerGroup; y++) {
                SAMRecord record = ArtificialSAMUtils.createArtificialRead(getHeader(), "readUno", 0, ++alignmentStart, 20);
                record.setAttribute("RG", groupRecord.getReadGroupId());
                records.add(record);
            }
        }

        List<String> filterList = new ArrayList<String>();
        filterList.add(privateTestDir + "readgroupblacklisttest.txt");

        ReadGroupBlackListFilter filter = new ReadGroupBlackListFilter(filterList);
        int filtered = 0;
        int unfiltered = 0;
        for (SAMRecord record : records) {
            String readGroup = record.getReadGroup().getReadGroupId();
            if (filter.filterOut(record)) {
                if (!("ReadGroup3".equals(readGroup) || "ReadGroup4".equals(readGroup)))
                    Assert.fail("Read group " + readGroup + " was filtered");
                filtered++;
            } else {
                if ("ReadGroup3".equals(readGroup) || "ReadGroup4".equals(readGroup))
                    Assert.fail("Read group " + readGroup + " was not filtered");
                unfiltered++;
            }
        }

        int filteredExpected = recordsPerGroup * 2;
        int unfilteredExpected = recordsPerGroup * (getReadGroupCount() - 2);
        Assert.assertEquals(filtered, filteredExpected, "Filtered");
        Assert.assertEquals(unfiltered, unfilteredExpected, "Uniltered");
    }

    @Test
    public void testFilterOutByListFile() {
        int recordsPerGroup = 3;
        List<SAMRecord> records = new ArrayList<SAMRecord>();
        int alignmentStart = 0;
        for (int x = 1; x <= getReadGroupCount(); x++) {
            SAMReadGroupRecord groupRecord = getHeader().getReadGroup(getReadGroupId(x));
            for (int y = 1; y <= recordsPerGroup; y++) {
                SAMRecord record = ArtificialSAMUtils.createArtificialRead(getHeader(), "readUno", 0, ++alignmentStart, 20);
                record.setAttribute("RG", groupRecord.getReadGroupId());
                records.add(record);
            }
        }

        List<String> filterList = new ArrayList<String>();
        filterList.add(privateTestDir + "readgroupblacklisttestlist.txt");

        ReadGroupBlackListFilter filter = new ReadGroupBlackListFilter(filterList);
        int filtered = 0;
        int unfiltered = 0;
        for (SAMRecord record : records) {
            String readGroup = record.getReadGroup().getReadGroupId();
            if (filter.filterOut(record)) {
                if (!("ReadGroup3".equals(readGroup) || "ReadGroup4".equals(readGroup)))
                    Assert.fail("Read group " + readGroup + " was filtered");
                filtered++;
            } else {
                if ("ReadGroup3".equals(readGroup) || "ReadGroup4".equals(readGroup))
                    Assert.fail("Read group " + readGroup + " was not filtered");
                unfiltered++;
            }
        }

        int filteredExpected = recordsPerGroup * 2;
        int unfilteredExpected = recordsPerGroup * (getReadGroupCount() - 2);
        Assert.assertEquals(filtered, filteredExpected, "Filtered");
        Assert.assertEquals(unfiltered, unfilteredExpected, "Uniltered");
    }
}
