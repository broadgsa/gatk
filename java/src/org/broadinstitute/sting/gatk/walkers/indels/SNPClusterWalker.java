package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.VariationRod;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.Variation;

@WalkerName("SNPClusters")
@Requires(value={DataSource.REFERENCE},referenceMetaData={@RMD(name="snps",type= VariationRod.class)})
public class SNPClusterWalker extends RefWalker<GenomeLoc, GenomeLoc> {
    @Argument(fullName="windowSize", shortName="window", doc="window size for calculating clusters", required=false)
    int windowSize = 10;

    public void initialize() {
        if ( windowSize < 1)
            throw new RuntimeException("Window Size must be a positive integer");
    }

    public GenomeLoc map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        Variation snp = (Variation)tracker.lookup("snps", null);
        return (snp != null && snp.isSNP()) ? context.getLocation() : null;
    }

    public void onTraversalDone(GenomeLoc sum) {
        if ( sum != null && sum.getStart() != sum.getStop() )
            out.println(sum);
    }

    public GenomeLoc reduceInit() {
        return null;
    }

    public GenomeLoc reduce(GenomeLoc value, GenomeLoc sum) {
        // ignore non-SNP variants
        if ( value == null )
            return sum;

        // if we have no previous SNPs start with the new location
        if ( sum == null )
            return value;

        // if we hit a new contig, emit and start with the new location
        if ( sum.getContigIndex() != value.getContigIndex() ) {
            if ( sum.getStart() != sum.getStop() )
                out.println(sum);
            return value;
        }

        // if the last SNP location was within a window, merge them
        if ( value.getStart() - sum.getStop() <= windowSize ) {
            sum = GenomeLocParser.setStop(sum,value.getStart());            
            return sum;
        }

        // otherwise, emit and start with the new location
        if ( sum.getStart() != sum.getStop() )
            out.println(sum);
        return value;
    }
}
