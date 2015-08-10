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

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.ArrayList;
import java.util.List;

/**
 * A class for analyzing and validating the read stream produced by an ArtificialSingleSampleReadStream.
 *
 * Collects various statistics about the stream of reads it's fed, and validates the stream
 * by checking whether the collected statistics match the nominal properties of the stream.
 *
 * Subclasses are expected to override the validate() method in order to check whether an artificial
 * read stream has been *transformed* in some way (eg., by downsampling or some other process), rather
 * than merely checking whether the stream matches its original properties.
 *
 * Usage is simple:
 *
 * ArtificialSingleSampleReadStreamAnalyzer analyzer = new ArtificialSingleSampleReadStreamAnalyzer(originalStream);
 * analyzer.analyze(originalOrTransformedStream);
 * analyzer.validate();  // override this method if you want to check whether the stream has been transformed
 *                       // in a certain way relative to the original stream
 *
 * @author David Roazen
 */
public class ArtificialSingleSampleReadStreamAnalyzer {
    protected ArtificialSingleSampleReadStream originalStream;
    protected SAMRecord lastRead;
    protected int totalReads;
    protected boolean allSamplesMatch;
    protected int numContigs;
    protected List<Integer> stacksPerContig;
    protected Integer minReadsPerStack;
    protected Integer maxReadsPerStack;
    protected Integer minDistanceBetweenStacks;
    protected Integer maxDistanceBetweenStacks;
    protected Integer minReadLength;
    protected Integer maxReadLength;
    protected int numUnmappedReads;

    protected int currentContigNumStacks;
    protected int currentStackNumReads;

    /**
     * Construct a new read stream analyzer, providing an ArtificialSingleSampleReadStream that will
     * serve as the basis for comparison after the analysis is complete.
     *
     * @param originalStream the original ArtificialSingleSampleReadStream upon which the stream
     *                       that will be fed to the analyzer is based
     */
    public ArtificialSingleSampleReadStreamAnalyzer( ArtificialSingleSampleReadStream originalStream ) {
        this.originalStream = originalStream;
        reset();
    }

    /**
     * Reset all read stream statistics collected by this analyzer to prepare for a fresh run
     */
    public void reset() {
        lastRead = null;
        totalReads = 0;
        allSamplesMatch = true;
        numContigs = 0;
        stacksPerContig = new ArrayList<Integer>();
        minReadsPerStack = null;
        maxReadsPerStack = null;
        minDistanceBetweenStacks = null;
        maxDistanceBetweenStacks = null;
        minReadLength = null;
        maxReadLength = null;
        numUnmappedReads = 0;
        currentContigNumStacks = 0;
        currentStackNumReads = 0;
    }

    /**
     * Collect statistics on the stream of reads passed in
     *
     * @param stream the stream of reads to analyze
     */
    public void analyze( Iterable<SAMRecord> stream ) {
        for ( SAMRecord read : stream ) {
            update(read);
        }
        finalizeStats();
    }

    /**
     * Validate the stream by checking whether our collected statistics match the properties of the
     * original stream. Throws a ReviewedGATKException if the stream is invalid.
     *
     * Override this method if you want to check whether the stream has been transformed in some
     * way relative to the original stream.
     */
    public void validate() {
        if ( (originalStream.getNumContigs() == 0 || originalStream.getNumStacksPerContig() == 0) && originalStream.getNumUnmappedReads() == 0 ) {
            if ( totalReads != 0 ) {
                throw new ReviewedGATKException("got reads from the stream, but the stream was configured to have 0 reads");
            }
            return;  // no further validation needed for the 0-reads case
        }
        else if ( totalReads == 0 ) {
            throw new ReviewedGATKException("got no reads from the stream, but the stream was configured to have > 0 reads");
        }

        if ( ! allSamplesMatch ) {
            throw new ReviewedGATKException("some reads had the wrong sample");
        }

        if ( numContigs != originalStream.getNumContigs() ) {
            throw new ReviewedGATKException("number of contigs not correct");
        }

        if ( stacksPerContig.size() != originalStream.getNumContigs() ) {
            throw new ReviewedGATKException(String.format("bug in analyzer code: calculated sizes for %d contigs even though there were only %d contigs",
                                                           stacksPerContig.size(), originalStream.getNumContigs()));
        }

        for ( int contigStackCount : stacksPerContig ) {
            if ( contigStackCount != originalStream.getNumStacksPerContig() ) {
                throw new ReviewedGATKException("contig had incorrect number of stacks");
            }
        }

        if ( originalStream.getNumStacksPerContig() > 0 ) {
            if ( minReadsPerStack < originalStream.getMinReadsPerStack() ) {
                throw new ReviewedGATKException("stack had fewer than the minimum number of reads");
            }
            if ( maxReadsPerStack > originalStream.getMaxReadsPerStack() ) {
                throw new ReviewedGATKException("stack had more than the maximum number of reads");
            }
        }
        else if ( minReadsPerStack != null || maxReadsPerStack != null ) {
            throw new ReviewedGATKException("bug in analyzer code: reads per stack was calculated even though 0 stacks per contig was specified");
        }

        if ( originalStream.getNumStacksPerContig() > 1 ) {
            if ( minDistanceBetweenStacks < originalStream.getMinDistanceBetweenStacks() ) {
                throw new ReviewedGATKException("stacks were separated by less than the minimum distance");
            }
            if ( maxDistanceBetweenStacks > originalStream.getMaxDistanceBetweenStacks() ) {
                throw new ReviewedGATKException("stacks were separated by more than the maximum distance");
            }
        }
        else if ( minDistanceBetweenStacks != null || maxDistanceBetweenStacks != null ) {
            throw new ReviewedGATKException("bug in analyzer code: distance between stacks was calculated even though numStacksPerContig was <= 1");
        }

        if ( minReadLength < originalStream.getMinReadLength() ) {
            throw new ReviewedGATKException("read was shorter than the minimum allowed length");
        }
        if ( maxReadLength > originalStream.getMaxReadLength() ) {
            throw new ReviewedGATKException("read was longer than the maximum allowed length");
        }

        if ( numUnmappedReads != originalStream.getNumUnmappedReads() ) {
            throw new ReviewedGATKException(String.format("wrong number of unmapped reads: requested %d but saw %d",
                                                           originalStream.getNumUnmappedReads(), numUnmappedReads));
        }

        if ( (originalStream.getNumContigs() == 0 || originalStream.getNumStacksPerContig() == 0) &&
             numUnmappedReads != totalReads ) {
            throw new ReviewedGATKException("stream should have consisted only of unmapped reads, but saw some mapped reads");
        }
    }

    public void update( SAMRecord read ) {
        if ( read.getReadUnmappedFlag() ) {
            numUnmappedReads++;

            if ( numUnmappedReads == 1 && lastRead != null ) {
                processContigChange();
                numContigs--;
            }
        }
        else if ( lastRead == null ) {
            numContigs = 1;
            currentContigNumStacks = 1;
            currentStackNumReads = 1;
        }
        else if ( ! read.getReferenceIndex().equals(lastRead.getReferenceIndex()) ) {
            processContigChange();
        }
        else if ( read.getAlignmentStart() != lastRead.getAlignmentStart() ) {
            processStackChangeWithinContig(read);
        }
        else {
            currentStackNumReads++;
        }

        updateReadLength(read.getReadLength());
        allSamplesMatch = allSamplesMatch && readHasCorrectSample(read);
        totalReads++;

        lastRead = read;
    }


    private void processContigChange() {
        numContigs++;

        stacksPerContig.add(currentContigNumStacks);
        currentContigNumStacks = 1;

        updateReadsPerStack(currentStackNumReads);
        currentStackNumReads = 1;
    }

    private void processStackChangeWithinContig( SAMRecord read ) {
        currentContigNumStacks++;

        updateReadsPerStack(currentStackNumReads);
        currentStackNumReads = 1;

        updateDistanceBetweenStacks(read.getAlignmentStart() - lastRead.getAlignmentStart());
    }

    private void updateReadsPerStack( int stackReadCount ) {
        if ( minReadsPerStack == null || stackReadCount < minReadsPerStack ) {
            minReadsPerStack = stackReadCount;
        }
        if ( maxReadsPerStack == null || stackReadCount > maxReadsPerStack ) {
            maxReadsPerStack = stackReadCount;
        }
    }

    private void updateDistanceBetweenStacks( int stackDistance ) {
        if ( minDistanceBetweenStacks == null || stackDistance < minDistanceBetweenStacks ) {
            minDistanceBetweenStacks = stackDistance;
        }
        if ( maxDistanceBetweenStacks == null || stackDistance > maxDistanceBetweenStacks ) {
            maxDistanceBetweenStacks = stackDistance;
        }
    }

    private void updateReadLength( int readLength ) {
        if ( minReadLength == null || readLength < minReadLength ) {
            minReadLength = readLength;
        }
        if ( maxReadLength == null || readLength > maxReadLength ) {
            maxReadLength = readLength;
        }
    }

    private boolean readHasCorrectSample( SAMRecord read ) {
        return originalStream.getReadGroupID().equals(read.getAttribute("RG"));
    }

    public void finalizeStats() {
        if ( lastRead != null && ! lastRead.getReadUnmappedFlag() ) {
            stacksPerContig.add(currentContigNumStacks);
            updateReadsPerStack(currentStackNumReads);
        }
    }
}
