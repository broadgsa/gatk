/*
* Copyright (c) 2012 The Broad Institute
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

import htsjdk.samtools.*;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.ReadProperties;
import org.broadinstitute.gatk.utils.ValidationExclusion;
import org.broadinstitute.gatk.engine.datasources.reads.SAMDataSource;
import org.broadinstitute.gatk.utils.exceptions.UserException;

/**
 * Filter out reads that are over-soft-clipped
 *
 * <p>
 *     This filter is intended to filter out reads that are potentially from foreign organisms.
 *     From experience with sequencing of human DNA we have found cases of contamination by bacterial
 *     organisms; the symptoms of such contamination are a class of reads with only a small number
 *     of aligned bases and additionally many soft-clipped bases on both ends.  This filter is intended
 *     to remove such reads.
 * </p>
 *
 */
public class OverclippedReadFilter extends ReadFilter {

    @Argument(fullName = "filter_is_too_short_value", shortName = "filterTooShort", doc = "Value for which reads with less than this number of aligned bases is considered too short", required = false)
    int tooShort = 30;


    public boolean filterOut(final SAMRecord read) {
        boolean sawLeadingSoftclip = false;
        boolean sawAlignedBase = false;
        int alignedLength = 0;

        for ( final CigarElement element : read.getCigar().getCigarElements() ) {
            if ( element.getOperator() == CigarOperator.S ) {
                if ( sawAlignedBase )  // if this is true then we must also have seen a leading soft-clip
                    return (alignedLength < tooShort);
                sawLeadingSoftclip = true;
            } else if ( element.getOperator().consumesReadBases() ) {   // M, I, X, and EQ (S was already accounted for above)
                if ( !sawLeadingSoftclip )
                    return false;
                sawAlignedBase = true;
                alignedLength += element.getLength();
            }
        }

        return false;
    }
}
