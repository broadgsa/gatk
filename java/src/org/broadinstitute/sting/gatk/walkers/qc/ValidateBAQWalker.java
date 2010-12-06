package org.broadinstitute.sting.gatk.walkers.qc;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.Argument;

import java.io.PrintStream;

/**
 * Walks over the input data set, calculating the number of reads seen for diagnostic purposes.
 * Can also count the number of reads matching a given criterion using read filters (see the
 * --read-filter command line argument).  Simplest example of a read-backed analysis.
 */
@BAQMode(QualityMode = BAQ.QualityMode.DONT_MODIFY, ApplicationTime = BAQ.ApplicationTime.HANDLED_IN_WALKER)
@Reference(window=@Window(start=-5,stop=5))
@Requires({DataSource.READS, DataSource.REFERENCE, DataSource.REFERENCE_BASES})
public class ValidateBAQWalker extends ReadWalker<Integer, Integer> {
    @Output(doc="File to which results should be written",required=true)
    protected PrintStream out;

    @Argument(doc="maximum read length to apply the BAQ calculation too",required=false)
    protected int maxReadLen = 1000;

    @Argument(doc="only operates on reads with this name",required=false)
    protected String readName = null;

    @Argument(doc="only prints out detailed information on the impact of BAQ when our implementation differences from the samtools BAQ tag", required=false)
    protected boolean onlyPrintOnFailures = false;

    @Argument(doc="Also prints out detailed comparison information when for known calculation differences", required=false)
    protected boolean alsoPrintWarnings = false;

    int counter = 0;

    BAQ baqHMM = new BAQ(1e-3, 0.1, 7, 0);         // matches current samtools parameters

    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker tracker) {
        IndexedFastaSequenceFile refReader = this.getToolkit().getReferenceDataSource().getReference();

        if ( (readName == null || readName.equals(read.getReadName())) && read.getReadLength() <= maxReadLen && BAQ.hasBAQTag(read) ) {
            byte[] baqFromTag = BAQ.calcBAQFromTag(read, false, false);
            if (counter++ % 1000 == 0) out.printf("Checking read %s (%d)%n", read.getReadName(), counter);
            BAQ.BAQCalculationResult baq = baqHMM.calcBAQFromHMM(read, refReader);

            boolean fail = false;
            boolean print = false;
            int badi = 0;
            for ( badi = 0; badi < baqFromTag.length; badi++ ) {
                if ( baqFromTag[badi] != baq.bq[badi] ) {
                    if (MathUtils.arrayMin(read.getBaseQualities()) == 0) {
                        print = true;
                        out.printf("  different, but Q0 base detected%n");
                        break;
                    }
                    else if (readHasSoftClip(read)) {
                        print = true;
                        out.printf("  different, but soft clip detected%n");
                        break;
                    } else if (readHasDeletion(read)) {
                        print = true;
                        out.printf("  different, but deletion detected%n");
                        break;
                    } else {
                        fail = true; print = true;
                        break;
                    }
                }
            }

            if ( ! onlyPrintOnFailures || ( print && ( alsoPrintWarnings || fail ) ) ) {
                out.printf("  read length   : %d%n", read.getReadLength());
                out.printf("  read start    : %d%n", read.getAlignmentStart());
                out.printf("  cigar         : %s%n", read.getCigarString());
                out.printf("  read bases    : %s%n", new String(read.getReadBases()));
                out.printf("  ref length    : %d%n", baq.refBases.length);
                out.printf("  BQ tag        : %s%n", read.getStringAttribute(BAQ.BAQ_TAG));
                printQuals("  BQ deltas     : ", getBAQDeltas(read), true);
                printQuals("  original quals: ", read.getBaseQualities(), true);
                printQuals("  original quals: ", read.getBaseQualities());
                printQuals("  tag      quals: ", baqFromTag);
                printQuals("  hmm      quals: ", baq.bq);
                out.printf("  read bases    : %s%n", new String(read.getReadBases()));
            }


            if ( fail )
                throw new StingException(String.format("BAQ from read and from HMM differ in read %s at position %d: tag qual = %d, hmm qual = %d",
                        read.getReadName(), badi, baqFromTag[badi], baq.bq[badi]));
        }

        return 1;
    }

    private final static boolean readHasSoftClip(SAMRecord read) {
        for (CigarElement e : read.getCigar().getCigarElements()) {
            if ( e.getOperator() == CigarOperator.SOFT_CLIP )
                return true;
        }

        return false;
    }

    private final static boolean readHasDeletion(SAMRecord read) {
        for (CigarElement e : read.getCigar().getCigarElements()) {
            if ( e.getOperator() == CigarOperator.DELETION )
                return true;
        }

        return false;
    }

    private final void printQuals( String prefix, byte[] quals ) {
        printQuals(prefix, quals, false);
    }

    private final void printQuals( String prefix, byte[] quals, boolean asInt ) {
        out.print(prefix);
        for ( int i = 0; i < quals.length; i++) {
            if ( asInt ) {
                out.print((int)quals[i]);
                if ( i+1 != quals.length ) out.print(",");
            } else
                out.print((char)(quals[i]+33));
        }
        out.println();
    }

    /**
     * Get the BAQ delta bytes from the tag in read.  Returns null if no BAQ tag is present.
     * @param read
     * @return
     */
    public static byte[] getBAQDeltas(SAMRecord read) {
        byte[] baq = BAQ.getBAQTag(read);
        if ( baq != null ) {
            byte[] deltas = new byte[baq.length];
            for ( int i = 0; i < deltas.length; i++)
                deltas[i] = (byte)(-1 * (baq[i] - 64));
            return deltas;
        } else
            return null;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}

