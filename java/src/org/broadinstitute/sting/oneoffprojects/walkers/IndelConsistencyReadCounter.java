/*
 * Copyright (c) 2010.
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

package org.broadinstitute.sting.oneoffprojects.walkers;

import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broad.tribble.util.variantcontext.VariantContext;

import java.util.*;

@By(DataSource.READS)
@Requires(value={},referenceMetaData=@RMD(name="indels", type=VariantContext.class))
// walker to count reads that are and are not consistent with homozygous indels
public class IndelConsistencyReadCounter extends ReadWalker<Integer, Integer> {

    private long consistentReads = 0, misalignedReads = 0;

    public boolean filter(ReferenceContext ref, SAMRecord read) {
        return !doNotTryToClean(read);
    }

    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        VariantContext indel = null;
        for ( Collection<GATKFeature> rods : metaDataTracker.getContigOffsetMapping().values() ) {
            Iterator<GATKFeature> rodIter = rods.iterator();
            while ( rodIter.hasNext() ) {
                Object rod = rodIter.next().getUnderlyingObject();
                if ( VariantContextAdaptors.canBeConvertedToVariantContext(rod)) {
                    VariantContext vc = VariantContextAdaptors.toVariantContext("", rod, ref);
                    if ( vc.getName().equals("indels") ) {
                        indel = vc;
                        break;
                    }
                }
            }
        }

        if ( indel != null ) {
            if ( read.getAlignmentEnd() == indel.getStart() )
                return 0;

            if ( !containsAnyIndel(read) || !containsIndel(read, indel) )
                misalignedReads++;
             else
                consistentReads++;
        }

        return 0;
    }

    private boolean doNotTryToClean(SAMRecord read) {
        return read.getReadUnmappedFlag() ||
                read.getNotPrimaryAlignmentFlag() ||
                read.getReadFailsVendorQualityCheckFlag() ||
                read.getMappingQuality() == 0 ||
                read.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START ||
                (BadMateFilter.hasBadMate(read));
    }

    private static boolean containsAnyIndel(final SAMRecord r) {
        final Cigar cigar = r.getCigar();
        if ( cigar == null )
            return false;

        for ( final CigarElement e : cigar.getCigarElements() ) {
        	if (e.getOperator() == CigarOperator.D || e.getOperator() == CigarOperator.I )
                return true;
        }

    	return false;
    }

    private static boolean containsIndel(final SAMRecord r, final VariantContext vc) {
        int indelStart = vc.getStart() + 1;
        int readPos = r.getAlignmentStart();

        if ( vc.isInsertion() && indelStart == readPos )
            return true;

        final Cigar cigar = r.getCigar();

        int idx = 0;
        while ( readPos < indelStart && idx < cigar.numCigarElements() ) {

            final CigarElement ce = cigar.getCigarElement(idx);
            switch ( ce.getOperator() ) {
                case M:
                case I:
                    readPos += ce.getLength();
                    break;
                default:
                    break;
            }

            idx++;
        }

        if ( idx == cigar.numCigarElements() )
            return false;

        if ( readPos != indelStart )
            return false;

        final CigarElement ce = cigar.getCigarElement(idx);
        if ( vc.isDeletion() )
            return ( ce.getOperator() == CigarOperator.D && ce.getLength() == vc.getReference().getBases().length);
        return ( ce.getOperator() == CigarOperator.I && ce.getLength() == vc.getAlternateAllele(0).getBases().length);
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {
        System.out.println(consistentReads + " reads were initially consistent");
        System.out.println(misalignedReads + " reads were initially misaligned");
    }
}