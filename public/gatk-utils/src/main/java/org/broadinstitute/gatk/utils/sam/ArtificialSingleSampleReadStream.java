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

package org.broadinstitute.gatk.utils.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.iterators.GATKSAMIterator;
import org.broadinstitute.gatk.utils.iterators.GATKSAMIteratorAdapter;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

/**
 * An artificial stream of reads from a single read group/sample with configurable characteristics
 * such as:
 *
 * -the number of contigs that the reads should be distributed across
 * -number of "stacks" of reads sharing the same alignment start position per contig
 * -the min/max number of reads in each stack (exact values chosen randomly from this range)
 * -the min/max distance between stack start positions (exact values chosen randomly from this range)
 * -the min/max length of each read (exact values chosen randomly from this range)
 * -the number of unmapped reads
 *
 * The cigar string for all reads will be *M, where * is the length of the read.
 *
 * @author David Roazen
 */
public class ArtificialSingleSampleReadStream implements Iterable<SAMRecord> {
    private SAMFileHeader header;
    private String readGroupID;
    private int numContigs;
    private int numStacksPerContig;
    private int minReadsPerStack;
    private int maxReadsPerStack;
    private int minDistanceBetweenStacks;
    private int maxDistanceBetweenStacks;
    private int minReadLength;
    private int maxReadLength;
    private int numUnmappedReads;

    private static final String READ_GROUP_TAG = "RG";

    public ArtificialSingleSampleReadStream( SAMFileHeader header,
                                             String readGroupID,
                                             int numContigs,
                                             int numStacksPerContig,
                                             int minReadsPerStack,
                                             int maxReadsPerStack,
                                             int minDistanceBetweenStacks,
                                             int maxDistanceBetweenStacks,
                                             int minReadLength,
                                             int maxReadLength,
                                             int numUnmappedReads ) {
        this.header = header;
        this.readGroupID = readGroupID;
        this.numContigs = numContigs;
        this.numStacksPerContig = numStacksPerContig;
        this.minReadsPerStack = minReadsPerStack;
        this.maxReadsPerStack = maxReadsPerStack;
        this.minDistanceBetweenStacks = minDistanceBetweenStacks;
        this.maxDistanceBetweenStacks = maxDistanceBetweenStacks;
        this.minReadLength = minReadLength;
        this.maxReadLength = maxReadLength;
        this.numUnmappedReads = numUnmappedReads;

        validateStreamParameters();
    }

    private void validateStreamParameters() {
        if ( header == null || readGroupID == null ) {
            throw new ReviewedGATKException("null SAMFileHeader or read group ID") ;
        }

        if ( header.getReadGroup(readGroupID) == null ) {
            throw new ReviewedGATKException(String.format("Read group %s not found in SAMFileHeader", readGroupID));
        }

        if ( numContigs < 0 || numStacksPerContig < 0 || minReadsPerStack < 0 || maxReadsPerStack < 0 ||
             minDistanceBetweenStacks < 0 || maxDistanceBetweenStacks < 0 || minReadLength < 0 || maxReadLength < 0 ||
             numUnmappedReads < 0 ) {
            throw new ReviewedGATKException("Read stream parameters must be >= 0");
        }

        if ( (numContigs == 0 && numStacksPerContig != 0) || (numContigs != 0 && numStacksPerContig == 0) ) {
            throw new ReviewedGATKException("numContigs and numStacksPerContig must either both be > 0, or both be 0");
        }

        if ( minReadsPerStack > maxReadsPerStack ) {
            throw new ReviewedGATKException("minReadsPerStack > maxReadsPerStack");
        }

        if ( minDistanceBetweenStacks > maxDistanceBetweenStacks ) {
            throw new ReviewedGATKException("minDistanceBetweenStacks > maxDistanceBetweenStacks");
        }

        if ( minReadLength > maxReadLength ) {
            throw new ReviewedGATKException("minReadLength > maxReadLength");
        }
    }

    public Iterator<SAMRecord> iterator() {
        return makeReads().iterator();
    }

    public GATKSAMIterator getGATKSAMIterator() {
        return GATKSAMIteratorAdapter.adapt(iterator());
    }

    public Collection<SAMRecord> makeReads() {
        Collection<SAMRecord> reads = new ArrayList<SAMRecord>(numContigs * numStacksPerContig * maxReadsPerStack);

        for ( int contig = 0; contig < numContigs; contig++ ) {
            int alignmentStart = 1;

            for ( int stack = 0; stack < numStacksPerContig; stack++ ) {
                reads.addAll(makeReadStack(contig, alignmentStart, MathUtils.randomIntegerInRange(minReadsPerStack, maxReadsPerStack)));
                alignmentStart += MathUtils.randomIntegerInRange(minDistanceBetweenStacks, maxDistanceBetweenStacks);
            }
        }

        if ( numUnmappedReads > 0 ) {
            reads.addAll(makeReadStack(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX, SAMRecord.NO_ALIGNMENT_START, numUnmappedReads));
        }

        return reads;
    }

    private Collection<SAMRecord> makeReadStack( int contig, int alignmentStart, int stackSize ) {
        Collection<SAMRecord> readStack = new ArrayList<SAMRecord>(stackSize);

        for ( int i = 0; i < stackSize; i++ ) {
            SAMRecord read = ArtificialSAMUtils.createArtificialRead(header,
                                                                     "foo",
                                                                     contig,
                                                                     alignmentStart,
                                                                     MathUtils.randomIntegerInRange(minReadLength, maxReadLength));
            read.setAttribute(READ_GROUP_TAG, readGroupID);
            readStack.add(read);
        }

        return readStack;
    }

    public SAMFileHeader getHeader() {
        return header;
    }

    public String getReadGroupID() {
        return readGroupID;
    }

    public int getNumContigs() {
        return numContigs;
    }

    public int getNumStacksPerContig() {
        return numStacksPerContig;
    }

    public int getMinReadsPerStack() {
        return minReadsPerStack;
    }

    public int getMaxReadsPerStack() {
        return maxReadsPerStack;
    }

    public int getMinDistanceBetweenStacks() {
        return minDistanceBetweenStacks;
    }

    public int getMaxDistanceBetweenStacks() {
        return maxDistanceBetweenStacks;
    }

    public int getMinReadLength() {
        return minReadLength;
    }

    public int getMaxReadLength() {
        return maxReadLength;
    }

    public int getNumUnmappedReads() {
        return numUnmappedReads;
    }
}
