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

package org.broadinstitute.gatk.engine.downsampling;

import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.sam.ArtificialSingleSampleReadStream;
import org.broadinstitute.gatk.utils.sam.ArtificialSingleSampleReadStreamAnalyzer;

/**
 * Class for analyzing an artificial read stream that has been positionally downsampled, and verifying
 * that the downsampling was done correctly without changing the stream in unexpected ways.
 *
 * @author David Roazen
 */
public class PositionallyDownsampledArtificialSingleSampleReadStreamAnalyzer extends ArtificialSingleSampleReadStreamAnalyzer {
    private int targetCoverage;

    public PositionallyDownsampledArtificialSingleSampleReadStreamAnalyzer( ArtificialSingleSampleReadStream originalStream, int targetCoverage ) {
        super(originalStream);
        this.targetCoverage = targetCoverage;
    }

    /**
     * Overridden validate() method that checks for the effects of positional downsampling in addition to checking
     * for whether the original properties of the stream not affected by downsampling have been preserved
     */
    @Override
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

            // Check for the effects of positional downsampling:
            int stackMinimumAfterDownsampling = Math.min(targetCoverage, originalStream.getMinReadsPerStack());
            int stackMaximumAfterDownsampling = targetCoverage;

            if ( minReadsPerStack < stackMinimumAfterDownsampling ) {
                throw new ReviewedGATKException("stack had fewer than the minimum number of reads after downsampling");
            }
            if ( maxReadsPerStack > stackMaximumAfterDownsampling ) {
                throw new ReviewedGATKException("stack had more than the maximum number of reads after downsampling");
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
}
