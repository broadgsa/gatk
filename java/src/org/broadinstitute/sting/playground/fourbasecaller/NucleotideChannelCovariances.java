package org.broadinstitute.sting.playground.fourbasecaller;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import org.broadinstitute.sting.playground.illumina.FourIntensity;

public class NucleotideChannelCovariances {
    private DoubleMatrix2D[] covs;
    private int[] counts;

    public NucleotideChannelCovariances() {
        counts = new int[4];

        covs = new DoubleMatrix2D[4];
        for (int i = 0; i < 4; i++) {
            covs[i] = (DoubleFactory2D.dense).make(4, 4);
        }
    }

    public void add(Nucleotide nuc, FourIntensity sig, NucleotideChannelMeans mus) {
        FourIntensity f = new FourIntensity(sig);
        f.subtract(mus.channelMeans(nuc));

        DoubleMatrix2D cov = covs[nuc.ordinal()];

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                cov.set(i, j, cov.get(i, j) + f.getChannelIntensity(i)*f.getChannelIntensity(j));
            }
        }
        
        counts[nuc.ordinal()]++;
    }

    public DoubleMatrix2D channelCovariances(Nucleotide base) {
        DoubleMatrix2D newcov = covs[base.ordinal()].copy();

        /*
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                newcov.set(i, j, newcov.get(i, j)/((double) counts[base.ordinal()]));
            }
        }
        */

        return newcov;
    }

    public void invert() {
        Algebra alg = new Algebra();

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    covs[i].set(j, k, covs[i].get(j, k)/((double) counts[i]));
                }
            }
            
            DoubleMatrix2D covinv = alg.inverse(covs[i]);
            covs[i] = covinv;
        }
    }
}
