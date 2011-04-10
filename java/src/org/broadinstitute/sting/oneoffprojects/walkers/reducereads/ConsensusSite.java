package org.broadinstitute.sting.oneoffprojects.walkers.reducereads;

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

/**
* Created by IntelliJ IDEA.
* User: depristo
* Date: 4/8/11
* Time: 3:01 PM
* To change this template use File | Settings | File Templates.
*/
final class ConsensusSite {
    final GenomeLoc loc;
    final Set<PileupElement> overlappingReads = new HashSet<PileupElement>();
    final int offset;
    final BaseCounts counts = new BaseCounts();

    ConsensusType markedType = null;

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

    public void setMarkedType(ConsensusType markedType) {
        this.markedType = markedType;
    }

    public ConsensusType getMarkedType() {
        if ( markedType == null ) throw new ReviewedStingException("markedType not yet set!");
        return markedType;
    }


}
