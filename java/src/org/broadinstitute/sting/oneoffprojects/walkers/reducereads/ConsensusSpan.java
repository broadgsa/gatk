package org.broadinstitute.sting.oneoffprojects.walkers.reducereads;

import org.broadinstitute.sting.utils.GenomeLoc;

/**
* Created by IntelliJ IDEA.
* User: depristo
* Date: 4/8/11
* Time: 3:01 PM
* To change this template use File | Settings | File Templates.
*/
final class ConsensusSpan {
    final int refStart; // the start position on the reference for relative calculations
    final GenomeLoc loc;
    final ConsensusType consensusType;

    public ConsensusSpan(final int refStart, GenomeLoc loc, ConsensusType consensusType) {
        this.refStart = refStart;
        this.loc = loc;
        this.consensusType = consensusType;
    }

    public int getOffsetFromStartOfSites() {
        return loc.getStart() - refStart;
    }

    public int getGenomeStart() {
        return loc.getStart();
    }

    public int getGenomeStop() {
        return loc.getStop();
    }

    public ConsensusType getConsensusType() {
        return consensusType;
    }

    public int size() {
        return getGenomeStop() - getGenomeStart() + 1;
    }

    public boolean isConserved() { return getConsensusType() == ConsensusType.CONSERVED; }
    public boolean isVariable() { return getConsensusType() == ConsensusType.VARIABLE; }

    public String toString() {
        return String.format("%s %s", consensusType, loc);
    }
}
