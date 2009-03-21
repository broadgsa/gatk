package org.broadinstitute.sting.projects.fourbasecaller;

import org.broadinstitute.sting.illumina.FourIntensity;

public class NucleotideChannelMeans {
    private FourIntensity[] sigs;
    private int[] counts;

    public NucleotideChannelMeans() {
        counts = new int[4];
        sigs = new FourIntensity[4];

        for (int base = 0; base < 4; base++) {
            sigs[base] = new FourIntensity();
        }
    }

    public void add(Nucleotide base, FourIntensity intensity) {
        sigs[base.ordinal()].add(intensity);
        counts[base.ordinal()]++;
    }

    public FourIntensity channelMeans(Nucleotide base) {
        FourIntensity meansig = new FourIntensity(sigs[base.ordinal()]);
        meansig.divide(counts[base.ordinal()]);

        return meansig;
    }

    public float channelMean(Nucleotide base, int channel) {
        FourIntensity meansig = channelMeans(base);

        return meansig.getChannelIntensity(channel);
    }

    public String toString() {
        return "a:{" + channelMean(Nucleotide.A, 0) + " "  + channelMean(Nucleotide.A, 1) + " " + channelMean(Nucleotide.A, 2) + " " + channelMean(Nucleotide.A, 3) + "} " +
               "c:{" + channelMean(Nucleotide.C, 0) + " "  + channelMean(Nucleotide.C, 1) + " " + channelMean(Nucleotide.C, 2) + " " + channelMean(Nucleotide.C, 3) + "} " +
               "g:{" + channelMean(Nucleotide.G, 0) + " "  + channelMean(Nucleotide.G, 1) + " " + channelMean(Nucleotide.G, 2) + " " + channelMean(Nucleotide.G, 3) + "} " +
               "t:{" + channelMean(Nucleotide.T, 0) + " "  + channelMean(Nucleotide.T, 1) + " " + channelMean(Nucleotide.T, 2) + " " + channelMean(Nucleotide.T, 3) + "}";
    }
}
