/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.features.samread.SAMReadFeature;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.commandline.Argument;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Arrays;

import net.sf.samtools.SAMRecord;

/**
 * Checks a given downsampled pileup against the full pileup to ensure that the downsampled pileup could
 * possibly be a valid version of the full pileup.
 *
 * @author mhanna
 * @version 0.1
 */
public class DownsamplingValidationWalker extends LocusWalker<Integer,Long> {
    @Argument(fullName="max_expected_number_of_reads",shortName="menr",doc="The expected number of reads chosed by the downsampler.  Fewer than this number might be added to a given alignment start, but more than this should never be.",required=true)
    private int maxExpectedNumberOfReads = 0;

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        ReadBackedPileup pileup = context.getBasePileup();
        Collection<Object> allFeatures = tracker.getReferenceMetaData("reads");

        Collection<SAMReadFeature> unsampledReadsStartingAtThisLocus = new ArrayList<SAMReadFeature>();
        for(Object featureCandidate: allFeatures) {
            if(featureCandidate instanceof SAMReadFeature) {
                SAMReadFeature feature = (SAMReadFeature)featureCandidate;
                if(feature.getReferenceName().equals(ref.getLocus().getContig()) && feature.getAlignmentStart() == ref.getLocus().getStart())
                    unsampledReadsStartingAtThisLocus.add(feature);
            }
        }
        Collection<SAMRecord> sampledReadsStartingAtThisLocus = new ArrayList<SAMRecord>();
        for(SAMRecord read: pileup.getReads()) {
            if(read.getReferenceName().equals(ref.getLocus().getContig()) && read.getAlignmentStart() == ref.getLocus().getStart())
                sampledReadsStartingAtThisLocus.add(read);
        }

        int matchingReadsFound = 0;
        if(unsampledReadsStartingAtThisLocus.isEmpty()) {
            if(!sampledReadsStartingAtThisLocus.isEmpty())
                throw new ReviewedStingException("Downsampler hallucinated a read starting at locus "+ref.getLocus());
        }
        else {
            boolean foundMatch = false;
            for(SAMReadFeature unsampledRead: unsampledReadsStartingAtThisLocus) {
                for(SAMRecord sampledRead: sampledReadsStartingAtThisLocus) {
                    if(unsampledRead.getReadName().equals(sampledRead.getReadName()) &&
                            Arrays.equals(unsampledRead.getReadBases(),sampledRead.getReadBases())) {
                        foundMatch = true;
                        matchingReadsFound++;
                    }
                }
            }

            if(!foundMatch)
                throw new ReviewedStingException("Downsampler failed to include any read starting at locus "+ref.getLocus());

            if(matchingReadsFound > maxExpectedNumberOfReads)
                throw new ReviewedStingException("Downsampler found too many reads starting at locus "+ref.getLocus());
        }

        return matchingReadsFound;
    }

    // Given result of map function
    public Long reduceInit() { return 0L; }
    public Long reduce(Integer value, Long sum) {
        return value + sum;
    }

    public Long treeReduce(Long lhs, Long rhs ) {
        return lhs+rhs;
    }
}
