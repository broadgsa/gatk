package org.broadinstitute.sting.gatk.filters;

import org.testng.Assert;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMReadGroupRecord;

import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

public class ReadGroupBlackListFilterUnitTest extends BaseTest {
    private static final int READ_GROUP_COUNT = 5;
    private static final String READ_GROUP_PREFIX = "ReadGroup";
    private static final String SAMPLE_NAME_PREFIX = "Sample";
    private static final String PLATFORM_PREFIX = "Platform";
    private static final String PLATFORM_UNIT_PREFIX = "Lane";
    private static SAMFileHeader header;

    @BeforeClass
    public void beforeClass() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);

        List<String> readGroupIDs = new ArrayList<String>();
        List<String> sampleNames = new ArrayList<String>();

        for (int i = 1; i <= READ_GROUP_COUNT; i++) {
            readGroupIDs.add(READ_GROUP_PREFIX + i);
            sampleNames.add(SAMPLE_NAME_PREFIX + i);
        }

        ArtificialSAMUtils.createEnumeratedReadGroups(header, readGroupIDs, sampleNames);

        for (int i = 1; i <= READ_GROUP_COUNT; i++) {
            SAMReadGroupRecord groupRecord = header.getReadGroup(READ_GROUP_PREFIX + i);
            groupRecord.setAttribute("PL", PLATFORM_PREFIX + (((i-1)%2)+1));
            groupRecord.setAttribute("PU", PLATFORM_UNIT_PREFIX + (((i-1)%3)+1));
        }
    }

    @Test(expectedExceptions=ReviewedStingException.class)
    public void testBadFilter() {
        List<String> badFilters = Collections.singletonList("bad");
        new ReadGroupBlackListFilter(badFilters);
    }
    @Test(expectedExceptions=ReviewedStingException.class)
    public void testBadFilterTag() {
        List<String> badFilters = Collections.singletonList("bad:filter");
        new ReadGroupBlackListFilter(badFilters);
    }

    @Test(expectedExceptions=ReviewedStingException.class)
    public void testBadFilterFile() {
        List<String> badFilters = Collections.singletonList("/foo/bar/rgbl.txt");
        new ReadGroupBlackListFilter(badFilters);
    }

    @Test
    public void testFilterReadGroup() {
        SAMRecord filteredRecord = ArtificialSAMUtils.createArtificialRead(header, "readUno", 0, 1, 20);
        filteredRecord.setAttribute("RG", READ_GROUP_PREFIX + "1");

        SAMRecord unfilteredRecord = ArtificialSAMUtils.createArtificialRead(header, "readDos", 0, 2, 20);
        unfilteredRecord.setAttribute("RG", READ_GROUP_PREFIX + "2");

        List<String> filterList = new ArrayList<String>();
        filterList.add("RG:" + READ_GROUP_PREFIX + "1");

        ReadGroupBlackListFilter filter = new ReadGroupBlackListFilter(filterList);
        Assert.assertTrue(filter.filterOut(filteredRecord));
        Assert.assertFalse(filter.filterOut(unfilteredRecord));
    }

    @Test
    public void testFilterPlatformUnit() {
        SAMRecord filteredRecord = ArtificialSAMUtils.createArtificialRead(header, "readUno", 0, 1, 20);
        filteredRecord.setAttribute("RG", READ_GROUP_PREFIX + "1");

        SAMRecord unfilteredRecord = ArtificialSAMUtils.createArtificialRead(header, "readDos", 0, 2, 20);
        unfilteredRecord.setAttribute("RG", READ_GROUP_PREFIX + "2");

        List<String> filterList = new ArrayList<String>();
        filterList.add("PU:" + PLATFORM_UNIT_PREFIX + "1");

        ReadGroupBlackListFilter filter = new ReadGroupBlackListFilter(filterList);
        Assert.assertTrue(filter.filterOut(filteredRecord));
        Assert.assertFalse(filter.filterOut(unfilteredRecord));
    }

    @Test
    public void testFilterOutByReadGroup() {
        int recordsPerGroup = 3;
        List<SAMRecord> records = new ArrayList<SAMRecord>();
        int alignmentStart = 0;
        for (int x = 1; x <= READ_GROUP_COUNT; x++) {
            SAMReadGroupRecord groupRecord = header.getReadGroup(READ_GROUP_PREFIX + x);
            for (int y = 1; y <= recordsPerGroup; y++) {
                SAMRecord record = ArtificialSAMUtils.createArtificialRead(header, "readUno", 0, ++alignmentStart, 20);
                record.setAttribute("RG", groupRecord.getReadGroupId());
                records.add(record);
            }
        }

        List<String> filterList = new ArrayList<String>();
        filterList.add("RG:" + READ_GROUP_PREFIX + "1");
        filterList.add("RG:" + READ_GROUP_PREFIX + "3");

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
        int unfilteredExpected = recordsPerGroup * (READ_GROUP_COUNT - 2);
        Assert.assertEquals(filtered, filteredExpected, "Filtered");
        Assert.assertEquals(unfiltered, unfilteredExpected, "Uniltered");
    }

    @Test
    public void testFilterOutByAttribute() {
        int recordsPerGroup = 3;
        List<SAMRecord> records = new ArrayList<SAMRecord>();
        int alignmentStart = 0;
        for (int x = 1; x <= READ_GROUP_COUNT; x++) {
            SAMReadGroupRecord groupRecord = header.getReadGroup(READ_GROUP_PREFIX + x);
            for (int y = 1; y <= recordsPerGroup; y++) {
                SAMRecord record = ArtificialSAMUtils.createArtificialRead(header, "readUno", 0, ++alignmentStart, 20);
                record.setAttribute("RG", groupRecord.getReadGroupId());
                records.add(record);
            }
        }

        List<String> filterList = new ArrayList<String>();
        filterList.add("PU:" + PLATFORM_UNIT_PREFIX + "1");

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
        for (int x = 1; x <= READ_GROUP_COUNT; x++) {
            SAMReadGroupRecord groupRecord = header.getReadGroup(READ_GROUP_PREFIX + x);
            for (int y = 1; y <= recordsPerGroup; y++) {
                SAMRecord record = ArtificialSAMUtils.createArtificialRead(header, "readUno", 0, ++alignmentStart, 20);
                record.setAttribute("RG", groupRecord.getReadGroupId());
                records.add(record);
            }
        }

        List<String> filterList = new ArrayList<String>();
        filterList.add(validationDataLocation + "readgroupblacklisttest.txt");

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
        int unfilteredExpected = recordsPerGroup * (READ_GROUP_COUNT - 2);
        Assert.assertEquals(filtered, filteredExpected, "Filtered");
        Assert.assertEquals(unfiltered, unfilteredExpected, "Uniltered");
    }

    @Test
    public void testFilterOutByListFile() {
        int recordsPerGroup = 3;
        List<SAMRecord> records = new ArrayList<SAMRecord>();
        int alignmentStart = 0;
        for (int x = 1; x <= READ_GROUP_COUNT; x++) {
            SAMReadGroupRecord groupRecord = header.getReadGroup(READ_GROUP_PREFIX + x);
            for (int y = 1; y <= recordsPerGroup; y++) {
                SAMRecord record = ArtificialSAMUtils.createArtificialRead(header, "readUno", 0, ++alignmentStart, 20);
                record.setAttribute("RG", groupRecord.getReadGroupId());
                records.add(record);
            }
        }

        List<String> filterList = new ArrayList<String>();
        filterList.add(validationDataLocation + "readgroupblacklisttestlist.txt");

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
        int unfilteredExpected = recordsPerGroup * (READ_GROUP_COUNT - 2);
        Assert.assertEquals(filtered, filteredExpected, "Filtered");
        Assert.assertEquals(unfiltered, unfilteredExpected, "Uniltered");
    }
}
