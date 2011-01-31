package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators;

/**
 * Created by IntelliJ IDEA. User: kiran Date: Nov 29, 2010 Time: 3:25:59 PM To change this template use File | Settings
 * | File Templates.
 */
class NewPhaseStats {
    public int neitherPhased;
    public int onlyCompPhased;
    public int onlyEvalPhased;
    public int phasesAgree;
    public int phasesDisagree;

    public NewPhaseStats() {
        this.neitherPhased = 0;
        this.onlyCompPhased = 0;
        this.onlyEvalPhased = 0;
        this.phasesAgree = 0;
        this.phasesDisagree = 0;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Neither phased: " + neitherPhased + "\tOnly Comp: " + onlyCompPhased + "\tOnly Eval: " + onlyEvalPhased + "\tSame phase: " + phasesAgree + "\tOpposite phase: " + phasesDisagree);
        return sb.toString();
    }

    public static String[] getFieldNamesArray() {
        return new String[]{"total", "neither", "only_comp", "only_eval", "both", "match", "switch", "switch_rate"};
    }

    public Object getField(int index) {
        switch (index) {
            case (0):
                return (neitherPhased + onlyCompPhased + onlyEvalPhased + phasesAgree + phasesDisagree);
            case (1):
                return neitherPhased;
            case (2):
                return onlyCompPhased;
            case (3):
                return onlyEvalPhased;
            case (4):
                return (phasesAgree + phasesDisagree);
            case (5):
                return phasesAgree;
            case (6):
                return phasesDisagree;
            case (7):
                return ((phasesDisagree == 0) ? 0 : ((double) phasesDisagree) / (phasesAgree + phasesDisagree));
            default:
                return -1;
        }
    }
}
