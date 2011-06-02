package org.broadinstitute.sting.playground.gatk.walkers.reducereads;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
* Created by IntelliJ IDEA.
* User: depristo
* Date: 4/8/11
* Time: 3:01 PM
*
*/
final class ConsensusSite {
    final GenomeLoc loc;
    final Set<PileupElement> overlappingReads = new HashSet<PileupElement>();
    final int offset;
    final BaseCounts counts = new BaseCounts();

    ConsensusSpan.Type markedType = null;

    public ConsensusSite(GenomeLoc loc, int offset) {
        this.loc = loc;

        this.offset = offset;

    }

    public GenomeLoc getLoc() {
        return loc;
    }

    public Set<PileupElement> getOverlappingReads() {
        return overlappingReads;
    }

    public void addOverlappingRead(PileupElement elt) {
        overlappingReads.add(elt);
        counts.incr(elt.getBase());
    }

    public boolean isStrongConsensus(final double maxFractionDisagreeingBases) {
        int mostCommon = counts.countOfMostCommonBase();
        int total = counts.totalCounts();
        double fractionCommonBase = (1.0 * mostCommon) / total;
        return (1 - fractionCommonBase) < maxFractionDisagreeingBases;
    }

    public final static class ConsensusBase {
        byte base, qual;

        public byte getBase() {
            return base;
        }

        public byte getQual() {
            return qual;
        }

        public ConsensusBase(byte base, byte qual) {
            this.base = base;
            this.qual = qual;
        }
    }

    public ConsensusBase getConsensus() {
        byte base = counts.baseWithMostCounts();
        int qual = 0;

        for ( PileupElement p : overlappingReads ) {
            if ( p.getBase() == base )
                qual++;
        }

        return new ConsensusBase(base, QualityUtils.boundQual(qual, (byte)64));
    }

    public String toString() {
        return counts.toString();
    }

    public void setMarkedType(ConsensusSpan.Type markedType) {
        this.markedType = markedType;
    }

    public ConsensusSpan.Type getMarkedType() {
        if ( markedType == null ) throw new ReviewedStingException("markedType not yet set!");
        return markedType;
    }
}
