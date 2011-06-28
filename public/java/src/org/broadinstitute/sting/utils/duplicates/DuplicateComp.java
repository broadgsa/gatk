package org.broadinstitute.sting.utils.duplicates;

public class DuplicateComp {
    public int getQLarger() {
        return qLarger;
    }

    public void setQLarger(int qLarger) {
        this.qLarger = qLarger;
    }

    public int getQSmaller() {
        return qSmaller;
    }

    public void setQSmaller(int qSmaller) {
        this.qSmaller = qSmaller;
    }

    public boolean isMismatchP() {
        return mismatchP;
    }

    public void setMismatchP(boolean mismatchP) {
        this.mismatchP = mismatchP;
    }

    private int qLarger;
    private int qSmaller;
    private boolean mismatchP;

    public DuplicateComp(int qLarger, int qSmaller, boolean misMatchP) {
        this.qLarger = qLarger;
        this.qSmaller = qSmaller;
        this.mismatchP = misMatchP;
    }

    public String toString() {
        return String.format("%d %d %b", qLarger, qSmaller, mismatchP);
    }
}