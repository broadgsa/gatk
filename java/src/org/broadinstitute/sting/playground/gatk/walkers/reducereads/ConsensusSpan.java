package org.broadinstitute.sting.playground.gatk.walkers.reducereads;

import com.google.java.contract.*;
import org.broadinstitute.sting.utils.GenomeLoc;

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
 * Represents a span of a consensus region (conserved, or variable) on the reference genome.  Supports
 * either absolute or relative (refStart) positioning of the span.
 */
final class ConsensusSpan {

    /**
     * The type of an span is either conserved (little variability within the span) or
     * variable (too many differences among the reads in the span to determine the exact
     * haplotype sequence).
     */
    public enum Type {
        CONSERVED, VARIABLE;

        public static Type otherType(Type t) {
            switch ( t ) {
                case CONSERVED: return VARIABLE;
                case VARIABLE: return CONSERVED;
            }
            return CONSERVED;
        }
    }


    final int refStart; // the start position on the reference for relative calculations
    final GenomeLoc loc;
    final Type consensusType;

    @Requires({"refStart >= 0", "loc != null", "consensusType != null"})
    @Ensures({"this.refStart == refStart", "this.loc.equals(loc)", "this.consensusType.equals(consensusType)"})
    public ConsensusSpan(final int refStart, GenomeLoc loc, ConsensusSpan.Type consensusType) {
        if ( refStart < 0 ) throw new RuntimeException("RefStart must be greater than 0: " + refStart);
        if ( loc == null ) throw new RuntimeException("Loc must not be null");
        if ( consensusType == null ) throw new RuntimeException("ConsensusType must not be null");

        this.refStart = refStart;
        this.loc = loc;
        this.consensusType = consensusType;
    }

    public int getOffsetFromStartOfSites() {
        return loc.getStart() - refStart;
    }

    @Ensures("result >= 0")
    public int getGenomeStart() {
        return loc.getStart();
    }

    @Ensures("result >= 0")
    public int getGenomeStop() {
        return loc.getStop();
    }

    public ConsensusSpan.Type getConsensusType() {
        return consensusType;
    }

    @Ensures("result >= 0")
    public int size() {
        return getGenomeStop() - getGenomeStart() + 1;
    }

    @Ensures("result == !isVariable()")
    public boolean isConserved() { return getConsensusType() == Type.CONSERVED; }

    @Ensures("result == !isConserved()")
    public boolean isVariable() { return getConsensusType() == Type.VARIABLE; }

    @Ensures("result != null")
    public String toString() {
        return String.format("%s %s", consensusType, loc);
    }
}
