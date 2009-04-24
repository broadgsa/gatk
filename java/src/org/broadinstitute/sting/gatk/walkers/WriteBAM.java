package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.QualityUtils;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMFileHeader;
import edu.mit.broad.picard.reference.ReferenceSequence;

import java.util.ArrayList;
import java.util.Random;
import java.io.File;

/**
 * ReadErrorRateWalker assesses the error rate per read position ('cycle') by comparing the
 * read to its home on the reference and noting the mismatch rate.  It ignores reads with
 * indels in them, treats high and low-quality references bases the same, and does not count
 * ambiguous bases as mismatches.  It's also thread-safe, so you can process a slew of reads
 * in short order.
 *
 * @author Kiran Garimella
 */
public class IOCrusherWalker extends ReadWalker<SAMRecord, ArrayList<SAMFileWriter>> {
    @Argument(shortName="nWaysOut",required=false,defaultValue="1")
    public int nWaysOut;

    @Argument(shortName="readScaling",required=false,defaultValue="1")
    public float readScaling;

    @Argument(shortName="outputBase", required=true)
    public String outputBase;

    public long nReadsRead = 0;
    public long nReadsWritten = 0;

    /**
     *
     */
    public SAMRecord map(LocusContext context, SAMRecord read) {
        nReadsRead++;
        return read;
    }

    /**
     * 
     */
    public ArrayList<SAMFileWriter> reduceInit() {
        SAMFileWriterFactory fact = new SAMFileWriterFactory();
        ArrayList<SAMFileWriter> outputs = new ArrayList<SAMFileWriter>(nWaysOut);
        for ( int i = 0; i < nWaysOut; i++ ) {
            SAMFileHeader header = this.getToolkit().getSamReader().getFileHeader();
            outputs.add(fact.makeBAMWriter(header, true, new File(outputBase + "." + i + ".bam")));
        }
        return outputs;
    }

    /**
     * Summarize the error rate data.
     *
     */
    public ArrayList<SAMFileWriter> reduce(SAMRecord read, ArrayList<SAMFileWriter> outputs) {
        for ( SAMFileWriter out : outputs ) {
            if ( readScaling >= 1.0 ) {
                int nCopies = (int)Math.ceil(readScaling);
                for ( int i = 0; i < nCopies; i++) {
                    out.addAlignment(read);
                    nReadsWritten++;
                }
            } else if ( Math.random() < readScaling ) {
                out.addAlignment(read);
                nReadsWritten++;
            }
        }
        
        return outputs;
    }

    /**
     *
     */
    public void onTraversalDone(ArrayList<SAMFileWriter> outputs) {
        for ( SAMFileWriter out : outputs ) {
            out.close();
        }
        System.out.printf("Reads: read %d written %d%n", nReadsRead, nReadsWritten);
    }
}