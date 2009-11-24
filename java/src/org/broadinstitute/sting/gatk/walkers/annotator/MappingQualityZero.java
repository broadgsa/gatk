package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.genotype.Genotype;
import net.sf.samtools.SAMRecord;

import java.util.List;


public class MappingQualityZero extends StandardVariantAnnotation {

    public Pair<String, String> annotate(ReferenceContext ref, ReadBackedPileup pileup, List<Genotype> genotypes) {
        List<SAMRecord> reads = pileup.getReads();
        int MQ0Count = 0;
        for (int i=0; i < reads.size(); i++) {
            if ( reads.get(i).getMappingQuality() == 0 )
                MQ0Count++;
        }
        return new Pair<String, String>("MAPQ0", String.format("%d", MQ0Count));
    }

    public boolean useZeroQualityReads() { return true; }
}