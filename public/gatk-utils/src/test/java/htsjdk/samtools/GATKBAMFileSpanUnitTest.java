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

package htsjdk.samtools;

import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Tests of functionality of union, intersection operators.
 */
public class GATKBAMFileSpanUnitTest {
    @Test
    public void testUnionOfEmptyFileSpans() {
        GATKBAMFileSpan empty1 = new GATKBAMFileSpan();
        GATKBAMFileSpan empty2 = new GATKBAMFileSpan();
        GATKBAMFileSpan union = empty1.union(empty2);
        Assert.assertEquals(union.getGATKChunks().size(),0,"Elements inserted in union of two empty sets");
    }

    @Test
    public void testUnionOfNonOverlappingFileSpans() {
        GATKBAMFileSpan regionOne = new GATKBAMFileSpan(new GATKChunk(0,65535));
        GATKBAMFileSpan regionTwo = new GATKBAMFileSpan(new GATKChunk(1<<16,(1<<16)|65535));
        GATKBAMFileSpan union = regionOne.union(regionTwo);
        Assert.assertEquals(union.getGATKChunks().size(),2,"Discontiguous elements were merged");
        Assert.assertEquals(union.getGATKChunks().get(0),regionOne.getGATKChunks().get(0),"Wrong chunk was first in list");
        Assert.assertEquals(union.getGATKChunks().get(1),regionTwo.getGATKChunks().get(0),"Wrong chunk was second in list");
    }

    @Test
    public void testUnionOfContiguousFileSpans() {
        // Region 1 ends at position adjacent to Region 2 start:
        // |---1----|---2----|

        GATKBAMFileSpan regionOne = new GATKBAMFileSpan(new GATKChunk(0,1<<16));
        GATKBAMFileSpan regionTwo = new GATKBAMFileSpan(new GATKChunk(1<<16,(1<<16)|65535));
        GATKBAMFileSpan union = regionOne.union(regionTwo);
        Assert.assertEquals(union.getGATKChunks().size(),1,"Elements to be merged were not.");
        Assert.assertEquals(union.getGATKChunks().get(0),new GATKChunk(0,(1<<16)|65535));
    }

    @Test
    public void testUnionOfFileSpansFirstRegionEndsWithinSecondRegion() {
        // Region 1 ends within Region 2:
        //        |---2----|
        // |---1----|

        GATKBAMFileSpan regionOne = new GATKBAMFileSpan(new GATKChunk(0,(1<<16)|32767));
        GATKBAMFileSpan regionTwo = new GATKBAMFileSpan(new GATKChunk(1<<16,(1<<16)|65535));
        GATKBAMFileSpan union = regionOne.union(regionTwo);
        Assert.assertEquals(union.getGATKChunks().size(),1,"Elements to be merged were not.");
        Assert.assertEquals(union.getGATKChunks().get(0),new GATKChunk(0,(1<<16)|65535));
    }

    @Test
    public void testUnionOfFileSpansFirstRegionEndsAtSecondRegionEnd() {
        // Region 1 ends at Region 2 end:
        //        |---2----|
        // |---1-----------|

        GATKBAMFileSpan regionOne = new GATKBAMFileSpan(new GATKChunk(0,(1<<16)|65535));
        GATKBAMFileSpan regionTwo = new GATKBAMFileSpan(new GATKChunk(1<<16,(1<<16)|65535));
        GATKBAMFileSpan union = regionOne.union(regionTwo);
        Assert.assertEquals(union.getGATKChunks().size(),1,"Elements to be merged were not.");
        Assert.assertEquals(union.getGATKChunks().get(0),new GATKChunk(0,(1<<16)|65535));
    }

    @Test
    public void testUnionOfFileSpansFirstRegionEndsAfterSecondRegionEnd() {
        // Region 1 ends after Region 2 end:
        //        |---2----|
        // |---1---------------|

        GATKBAMFileSpan regionOne = new GATKBAMFileSpan(new GATKChunk(0,(1<<16)|65535));
        GATKBAMFileSpan regionTwo = new GATKBAMFileSpan(new GATKChunk(1<<16,(1<<16)|32767));
        GATKBAMFileSpan union = regionOne.union(regionTwo);
        Assert.assertEquals(union.getGATKChunks().size(),1,"Elements to be merged were not.");
        Assert.assertEquals(union.getGATKChunks().get(0),new GATKChunk(0,(1<<16)|65535));
    }

    @Test
    public void testUnionOfFileSpansFirstRegionStartsAtSecondRegionStart() {
        // Region 1 starts at Region 2 start, but ends before Region 2:
        // |---2--------|
        // |---1----|

        GATKBAMFileSpan regionOne = new GATKBAMFileSpan(new GATKChunk(1<<16,(1<<16)|32767));
        GATKBAMFileSpan regionTwo = new GATKBAMFileSpan(new GATKChunk(1<<16,(1<<16)|65535));
        GATKBAMFileSpan union = regionOne.union(regionTwo);
        Assert.assertEquals(union.getGATKChunks().size(),1,"Elements to be merged were not.");
        Assert.assertEquals(union.getGATKChunks().get(0),new GATKChunk(1<<16,(1<<16)|65535));
    }

    @Test
    public void testUnionOfFileSpansFirstRegionEqualToSecondRegion() {
        // Region 1 and Region 2 represent the same region:
        // |---2----|
        // |---1----|

        GATKBAMFileSpan regionOne = new GATKBAMFileSpan(new GATKChunk(1<<16,(1<<16)|65535));
        GATKBAMFileSpan regionTwo = new GATKBAMFileSpan(new GATKChunk(1<<16,(1<<16)|65535));
        GATKBAMFileSpan union = regionOne.union(regionTwo);
        Assert.assertEquals(union.getGATKChunks().size(),1,"Elements to be merged were not.");
        Assert.assertEquals(union.getGATKChunks().get(0),new GATKChunk(1<<16,(1<<16)|65535));
    }

    @Test
    public void testUnionOfStringOfFileSpans() {
        GATKBAMFileSpan regionOne = new GATKBAMFileSpan(new GATKChunk[] { new GATKChunk(0,1<<16), new GATKChunk(2<<16,3<<16) });
        GATKBAMFileSpan regionTwo = new GATKBAMFileSpan(new GATKChunk(1<<16,2<<16));
        GATKBAMFileSpan union = regionOne.union(regionTwo);
        Assert.assertEquals(union.getGATKChunks().size(),1,"Elements to be merged were not.");
        Assert.assertEquals(union.getGATKChunks().get(0),new GATKChunk(0,3<<16));
    }

    @Test
    public void testUnionAllFileSpansAdded() {
        GATKBAMFileSpan regionOne = new GATKBAMFileSpan(new GATKChunk[] { new GATKChunk(0,1<<16), new GATKChunk(2<<16,3<<16), new GATKChunk(20<<16,21<<16) });
        GATKBAMFileSpan regionTwo = new GATKBAMFileSpan(new GATKChunk(1<<16,2<<16));
        GATKBAMFileSpan union = regionOne.union(regionTwo);
        Assert.assertEquals(union.getGATKChunks().size(),2,"Elements to be merged were not.");
        Assert.assertEquals(union.getGATKChunks().get(0),new GATKChunk(0,3<<16));
        Assert.assertEquals(union.getGATKChunks().get(1),new GATKChunk(20<<16,21<<16));
    }

    @Test
    public void testIntersectionOfEmptyFileSpans() {
        GATKBAMFileSpan empty1 = new GATKBAMFileSpan();
        GATKBAMFileSpan empty2 = new GATKBAMFileSpan();
        GATKBAMFileSpan intersection = empty1.intersection(empty2);
        Assert.assertEquals(intersection.getGATKChunks().size(),0,"Elements inserted in intersection of two empty sets");
    }

    @Test
    public void testIntersectionOfNonOverlappingFileSpans() {
        GATKBAMFileSpan regionOne = new GATKBAMFileSpan(new GATKChunk(0,1<<16));
        GATKBAMFileSpan regionTwo = new GATKBAMFileSpan(new GATKChunk(1<<16,2<<16));
        GATKBAMFileSpan intersection = regionOne.intersection(regionTwo);
        Assert.assertEquals(intersection.getGATKChunks().size(),0,"Elements inserted in intersection of two non-intersecting filespans");
    }

    @Test
    public void testIntersectionOfSmallOverlapInFileSpans() {
        GATKBAMFileSpan regionOne = new GATKBAMFileSpan(new GATKChunk(0,1<<16));
        GATKBAMFileSpan regionTwo = new GATKBAMFileSpan(new GATKChunk(65535,2<<16));
        GATKBAMFileSpan intersection = regionOne.intersection(regionTwo);
        Assert.assertEquals(intersection.getGATKChunks().size(),1,"No intersection found between two partially overlapping filespans");
        Assert.assertEquals(intersection.getGATKChunks().get(0),new GATKChunk(65535,1<<16),"Determined intersection is incorrect.");
    }

    @Test
    public void testIntersectionOfStrictSubset() {
        GATKBAMFileSpan regionOne = new GATKBAMFileSpan(new GATKChunk(0,1<<16));
        GATKBAMFileSpan regionTwo = new GATKBAMFileSpan(new GATKChunk(0,2<<16));
        GATKBAMFileSpan intersection = regionOne.intersection(regionTwo);
        Assert.assertEquals(intersection.getGATKChunks().size(),1,"No intersection found between two partially overlapping filespans");
        Assert.assertEquals(intersection.getGATKChunks().get(0),new GATKChunk(0<<16,1<<16),"Determined intersection is incorrect.");

        // Make sure intersection is symmetric
        intersection = regionTwo.intersection(regionOne);
        Assert.assertEquals(intersection.getGATKChunks().size(),1,"No intersection found between two partially overlapping filespans");
        Assert.assertEquals(intersection.getGATKChunks().get(0),new GATKChunk(0<<16,1<<16),"Determined intersection is incorrect.");
    }

    @Test
    public void testIntersectionOfPartialOverlap() {
        GATKBAMFileSpan regionOne = new GATKBAMFileSpan(new GATKChunk(0,2<<16));
        GATKBAMFileSpan regionTwo = new GATKBAMFileSpan(new GATKChunk(0<<16|32768,1<<16|32768));
        GATKBAMFileSpan intersection = regionOne.intersection(regionTwo);
        Assert.assertEquals(intersection.getGATKChunks().size(),1,"No intersection found between two partially overlapping filespans");
        Assert.assertEquals(intersection.getGATKChunks().get(0),new GATKChunk(0<<16|32768,1<<16|32768),"Determined intersection is incorrect.");
    }

    @Test
    public void testIntersectionOfChunkLists() {
        GATKBAMFileSpan regionOne = new GATKBAMFileSpan(new GATKChunk(0,5<<16));
        GATKBAMFileSpan regionTwo = new GATKBAMFileSpan(new GATKChunk[] { new GATKChunk(1<<16,2<<16), new GATKChunk(3<<16,4<<16) });
        GATKBAMFileSpan intersection = regionOne.intersection(regionTwo);
        Assert.assertEquals(intersection.getGATKChunks().size(),2,"Wrong number of intersections found.");
        Assert.assertEquals(intersection.getGATKChunks().get(0),new GATKChunk(1<<16,2<<16),"Determined intersection is incorrect.");
        Assert.assertEquals(intersection.getGATKChunks().get(1),new GATKChunk(3<<16,4<<16),"Determined intersection is incorrect.");

        // Make sure intersection is symmetric
        intersection = regionTwo.intersection(regionOne);
        Assert.assertEquals(intersection.getGATKChunks().size(),2,"Wrong number of intersections found.");
        Assert.assertEquals(intersection.getGATKChunks().get(0),new GATKChunk(1<<16,2<<16),"Determined intersection is incorrect.");
        Assert.assertEquals(intersection.getGATKChunks().get(1),new GATKChunk(3<<16,4<<16),"Determined intersection is incorrect.");
    }

    @Test
    public void testSubtractionOfEmptyChunkLists() {
        GATKBAMFileSpan regionOne = new GATKBAMFileSpan();
        GATKBAMFileSpan regionTwo = new GATKBAMFileSpan();
        GATKBAMFileSpan subtraction = regionOne.minus(regionTwo);
        Assert.assertEquals(subtraction.getGATKChunks().size(),0,"Elements inserted in subtraction of two empty sets");
    }

    @Test
    public void testSingleIntervalSubtractedAway() {
        GATKBAMFileSpan regionOne = new GATKBAMFileSpan(new GATKChunk(0,1<<16));
        GATKBAMFileSpan regionTwo = new GATKBAMFileSpan(new GATKChunk(0,1<<16));
        GATKBAMFileSpan subtraction = regionOne.minus(regionTwo);
        Assert.assertEquals(subtraction.getGATKChunks().size(),0,"Elements inserted in complete subtraction of region");
    }

    @Test
    public void testMultipleIntervalsSubtractedAway() {
        GATKBAMFileSpan regionOne = new GATKBAMFileSpan(new GATKChunk[] { new GATKChunk(0,1<<16), new GATKChunk(2<<16,3<<16) });
        GATKBAMFileSpan regionTwo = new GATKBAMFileSpan(new GATKChunk[] { new GATKChunk(0,1<<16), new GATKChunk(2<<16,3<<16) });
        GATKBAMFileSpan subtraction = regionOne.minus(regionTwo);
        Assert.assertEquals(subtraction.getGATKChunks().size(),0,"Elements inserted in complete subtraction of region");
    }

    @Test
    public void testSubtractionOfStrictSubset() {
        GATKBAMFileSpan regionOne = new GATKBAMFileSpan(new GATKChunk(0,2<<16));
        GATKBAMFileSpan regionTwo = new GATKBAMFileSpan(new GATKChunk(0,1<<16));
        GATKBAMFileSpan subtraction = regionOne.minus(regionTwo);
        Assert.assertEquals(subtraction.getGATKChunks().size(),1,"Incorrect size in strict subset subtraction of region");
        Assert.assertEquals(subtraction.getGATKChunks().get(0),new GATKChunk(1<<16,2<<16),"Determined subtraction is incorrect.");
    }

    @Test
    public void testSubtractionOfPartialOverlap() {
        GATKBAMFileSpan regionOne = new GATKBAMFileSpan(new GATKChunk(0,2<<16));
        GATKBAMFileSpan regionTwo = new GATKBAMFileSpan(new GATKChunk(1<<16,3<<16));
        GATKBAMFileSpan subtraction = regionOne.minus(regionTwo);
        Assert.assertEquals(subtraction.getGATKChunks().size(),1,"Incorrect size in partial subset subtraction of region");
        Assert.assertEquals(subtraction.getGATKChunks().get(0),new GATKChunk(0<<16,1<<16),"Determined subtraction is incorrect.");
    }
}
