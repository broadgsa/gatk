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
public abstract class LocalAssemblyEngine {

    public enum ASSEMBLER {
        SIMPLE_DE_BRUIJN
    }

    private PrintStream out;
    private IndexedFastaSequenceFile referenceReader;

    protected LocalAssemblyEngine(PrintStream out, IndexedFastaSequenceFile referenceReader) {
        this.out = out;
        this.referenceReader = referenceReader;
    }

    protected PrintStream getOutputStream() { return out; }

    protected IndexedFastaSequenceFile getReferenceReader() { return referenceReader; }

    public abstract void runLocalAssembly(List<SAMRecord> reads);

}
