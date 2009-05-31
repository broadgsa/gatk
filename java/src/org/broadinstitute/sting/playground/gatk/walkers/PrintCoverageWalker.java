package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.StringUtil;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.playground.indels.SWPairwiseAlignment;
import org.broadinstitute.sting.playground.utils.GenotypeLikelihoods;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

@By(DataSource.REFERENCE)
public class PrintCoverageWalker extends LocusWalker<Integer, Integer> {
    public void initialize() {
        out.println("chrom position reference a c g t");
    }

    public String walkerType() { return "ByLocus"; }

    // Do we actually want to operate on the context?
    public boolean filter(RefMetaDataTracker tracker, char ref, LocusContext context) {
        return true;    // We are keeping all the reads
    }

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        int aCount = 0;
        int cCount = 0;
        int gCount = 0;
        int tCount = 0;

        char upRef = Character.toUpperCase(ref);

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
            char base = read.getReadString().charAt(offset);
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