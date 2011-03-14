package org.broadinstitute.sting.playground.gatk.walkers.assembly;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMRecord;

import java.io.PrintStream;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Mar 14, 2011
 */
public class SimpleDeBruijnAssembler extends LocalAssemblyEngine {

    // minimum base quality required in a contiguous stretch of a given read to be used in the assembly
    private static final int MIN_BASE_QUAL_TO_USE = 20;

    // reference base padding size
    private static final int REFERENCE_PADDING = 30;

    public SimpleDeBruijnAssembler(PrintStream out, IndexedFastaSequenceFile referenceReader) {
        super(out, referenceReader);
    }

    public void runLocalAssembly(List<SAMRecord> reads) {
        createDeBruijnGraph(reads);
        
        assignReadsToGraph(reads);
    }

    private void createDeBruijnGraph(List<SAMRecord> reads) {

        // TODO -- implement me

    }

    private void assignReadsToGraph(List<SAMRecord> reads) {

        // TODO -- implement me

    }
}
