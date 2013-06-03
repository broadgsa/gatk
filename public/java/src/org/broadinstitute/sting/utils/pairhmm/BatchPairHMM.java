package org.broadinstitute.sting.utils.pairhmm;

import java.util.List;

import org.broadinstitute.sting.utils.haplotype.Haplotype;

public interface BatchPairHMM {
    public void batchAdd(final List<Haplotype> haplotypes, 
			 final byte[] readBases,
			 final byte[] readQuals,
			 final byte[] insertionGOP,
			 final byte[] deletionGOP,
			 final byte[] overallGCP);

    public double[] batchGetResult();
}
