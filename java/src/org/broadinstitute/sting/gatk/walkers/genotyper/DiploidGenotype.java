package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.Utils;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Aug 4, 2009
 * Time: 6:46:09 PM
 * To change this template use File | Settings | File Templates.
 */
public enum DiploidGenotype {
    AA,
    AC,
    AG,
    AT,
    CC,
    CG,
    CT,
    GG,
    GT,
    TT;

    public boolean isHomRef(char r) {
        return isHom() && r == this.toString().charAt(0);
    }

    public boolean isHomVar(char r) {
        return isHom() && r != this.toString().charAt(0);
    }

    public boolean isHetRef(char r) {
        return Utils.countOccurrences(r, this.toString()) == 1;
    }

    public boolean isHom() {
        return ! isHet();
    }

    public boolean isHet() {
        switch (this) {
            case AA:
            case CC:
            case GG:
            case TT: return false;
            default: return true;
        }
    }
}