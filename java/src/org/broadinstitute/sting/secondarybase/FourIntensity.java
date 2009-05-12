package org.broadinstitute.sting.secondarybase;

import java.util.StringTokenizer;

public class FourIntensity {
    private float[] fIntensities;

    public FourIntensity() {
        fIntensities = new float[4];
    }

    public FourIntensity(float[] fIntensities) {
        this.fIntensities = fIntensities;
    }

    public FourIntensity(FourIntensity intensity) {
        fIntensities = new float[4];
        
        for (int channel = 0; channel < 4; channel++) {
            fIntensities[channel] = intensity.getChannelIntensity(channel);
        }
    }

    public void add(FourIntensity intensity) {
        for (int channel = 0; channel < 4; channel++) {
            fIntensities[channel] += intensity.getChannelIntensity(channel);
        }
    }

    public void subtract(FourIntensity intensity) {
        for (int channel = 0; channel < 4; channel++) {
            fIntensities[channel] -= intensity.getChannelIntensity(channel);
        }
    }

    public void divide(float divisor) {
        for (int channel = 0; channel < 4; channel++) {
            fIntensities[channel] /= divisor;
        }
    }

    public float getChannelIntensity(int channel) { return fIntensities[channel]; }

    public int brightestChannel() {
        int brightest = 0;

        for (int channel = 1; channel < 4; channel++) {
           if (fIntensities[channel] > fIntensities[brightest]) {
               brightest = channel;
           }
        }

        return brightest;
    }

    public String toString() {
        return "("  + getChannelIntensity(0) +
               ", " + getChannelIntensity(1) +
               ", " + getChannelIntensity(2) +
               ", " + getChannelIntensity(3) +
               ")";
    }
}
