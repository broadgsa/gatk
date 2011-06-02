package org.broadinstitute.sting.playground.gatk.walkers.reducereads;

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
* To change this template use File | Settings | File Templates.
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

    public ConsensusSpan(final int refStart, GenomeLoc loc, ConsensusSpan.Type consensusType) {
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

    public ConsensusSpan.Type getConsensusType() {
        return consensusType;
    }

    public int size() {
        return getGenomeStop() - getGenomeStart() + 1;
    }

    public boolean isConserved() { return getConsensusType() == Type.CONSERVED; }
    public boolean isVariable() { return getConsensusType() == Type.VARIABLE; }

    public String toString() {
        return String.format("%s %s", consensusType, loc);
    }
}
