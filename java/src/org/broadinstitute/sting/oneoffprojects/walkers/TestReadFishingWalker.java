package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.alignment.bwa.BWAAligner;
import org.broadinstitute.sting.alignment.bwa.BWAConfiguration;
import org.broadinstitute.sting.alignment.bwa.c.BWACAligner;
import org.broadinstitute.sting.alignment.Alignment;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.StingException;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.picard.reference.ReferenceSequence;

import java.io.IOException;

/**
 * A walker to experiment with fishing for reads in the GATK.  Has very limited utility in its current state.
 *
 * @author mhanna
 * @version 0.1
 */
public class TestReadFishingWalker extends ReadWalker<Integer,Long> {
    /**
     * An aligner for the small custom reference.
     */
    private BWAAligner aligner;

    @Override
    public void initialize()  {
        ReferenceSequence initialBases;

        try {
            IndexedFastaSequenceFile reader = new IndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
            SAMSequenceRecord firstSequence = reader.getSequenceDictionary().getSequences().get(0);
            initialBases = reader.getSubsequenceAt(firstSequence.getSequenceName(),1,100);
        }
        catch(IOException ex) {
            throw new StingException("Unable to load initial bases from reference sequence");
        }

        aligner = new BWACAligner(initialBases.getBases(),new BWAConfiguration());
    }

    @Override
    public Integer map(char[] ref, SAMRecord read) {
        Alignment bestAlignment = aligner.getBestAlignment(read.getReadBases());
        System.out.println("bestAlignment = " + bestAlignment);
        return 1;
    }


    /**
     * Provide an initial value for reduce computations.
     * @return Initial value of reduce.
     */
    @Override
    public Long reduceInit() {
        return 0L;
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     * @param value result of the map.
     * @param accum accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    @Override
    public Long reduce(Integer value, Long accum) {
        return value + accum;
    }

    @Override
    public void onTraversalDone(Long result) {
        aligner.close();
        super.onTraversalDone(result);
    }

}
