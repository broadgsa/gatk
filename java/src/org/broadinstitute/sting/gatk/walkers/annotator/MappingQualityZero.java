package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.Variation;
import net.sf.samtools.SAMRecord;

import java.util.List;


public class MappingQualityZero extends StandardVariantAnnotation {

    public String annotate(ReferenceContext ref, ReadBackedPileup pileup, Variation variation, List<Genotype> genotypes) {
        List<SAMRecord> reads = pileup.getReads();
        int MQ0Count = 0;
        for (int i=0; i < reads.size(); i++) {
            if ( reads.get(i).getMappingQuality() == 0 )
                MQ0Count++;
        }
        return String.format("%d", MQ0Count);
    }

    public String getKeyName() { return "MQ0"; }

    public String getDescription() { return "MQ0,1,Integer,\"Total Mapping Quality Zero Reads\""; }

    public boolean useZeroQualityReads() { return true; }
}