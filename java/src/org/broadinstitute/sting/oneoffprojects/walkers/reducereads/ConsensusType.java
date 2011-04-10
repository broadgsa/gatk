package org.broadinstitute.sting.oneoffprojects.walkers.reducereads;

/**
* Created by IntelliJ IDEA.
* User: depristo
* Date: 4/9/11
* Time: 7:52 PM
* To change this template use File | Settings | File Templates.
*/
public enum ConsensusType {
    CONSERVED, VARIABLE;

    public static ConsensusType otherType(ConsensusType t) {
        switch ( t ) {
            case CONSERVED: return VARIABLE;
            case VARIABLE: return CONSERVED;
        }
        return CONSERVED;
    }
}
