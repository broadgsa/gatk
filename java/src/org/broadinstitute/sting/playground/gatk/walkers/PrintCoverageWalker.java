package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import java.util.List;

@By(DataSource.REFERENCE)
public class PrintCoverageWalker extends LocusWalker<Integer, Integer> {
    public void initialize() {
        out.println("chrom position reference a c g t");
    }

    public String walkerType() { return "ByLocus"; }

    // Do we actually want to operate on the context?
    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return true;    // We are keeping all the reads
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        int aCount = 0;
        int cCount = 0;
        int gCount = 0;
        int tCount = 0;

        char upRef = Character.toUpperCase(ref.getBase());

        List<SAMRecord> reads = context.getReads();
        for ( int i = 0; i < reads.size(); i++ )
        {
            SAMRecord read = reads.get(i);

            if (read.getNotPrimaryAlignmentFlag() ||
                read.getDuplicateReadFlag() ||
                read.getReadUnmappedFlag() ||
                read.getMappingQuality() <= 0
                    ) {
                continue;
            }

            int offset = context.getOffsets().get(i);
            char base = (char)read.getReadBases()[offset];
            if (base == 'a' || base == 'A') { aCount++; }
            if (base == 'c' || base == 'C') { cCount++; }
            if (base == 'g' || base == 'G') { gCount++; }
            if (base == 't' || base == 'T') { tCount++; }
        }


        out.printf("%s %d %s %d %d %d %d\n",
                   context.getContig(),
                   context.getPosition(),
                   upRef,
                   aCount,
                   cCount,
                   gCount,
                   tCount);

        return -1;
    }


    // Given result of map function
    public Integer reduceInit() {
        return 0;
    }
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    @Override
    public void onTraversalDone(Integer result) {
    }


}
