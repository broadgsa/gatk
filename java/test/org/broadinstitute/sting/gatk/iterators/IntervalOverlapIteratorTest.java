package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.junit.Assert;
import static org.junit.Assert.fail;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Iterator;


/**
 * @author aaron
 *         <p/>
 *         Class DuplicateDetectorIteratorTest
 *         <p/>
 *         test the DuplicateDetectorIterator class.
 */
public class IntervalOverlapIteratorTest extends BaseTest {
    private final File bam = new File(validationDataLocation + "index_test.bam");
    private static IndexedFastaSequenceFile seq;
    private int chromosomeOneReadCount = 885;

    //GenomeLoc.setupRefContigOrdering(seq.getSequenceDictionary());
    /**
     * This function does the setup of our parser, before each method call.
     * <p/>
     * Called before every test case method.
     */
    @BeforeClass
    public static void beforeAll() {
        try {
            seq = new IndexedFastaSequenceFile(new File(seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"));
        } catch (FileNotFoundException e) {
            fail("Unexpected Exception" + e.getLocalizedMessage());
        }
        GenomeLocParser.setupRefContigOrdering(seq);
    }

    @Test
    public void testOverlappingIntervals() {
        int countOfReads = 0;
        int seqLength = seq.getSequenceDictionary().getSequence("chr1").getSequenceLength();
        GenomeLoc last = GenomeLocParser.createGenomeLoc("chr1", 1, 470535);

        // first count the initial pile of reads
        SAMFileReader reader = new SAMFileReader(bam, true);
        reader.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
        Iterator<SAMRecord> i = reader.queryOverlapping("chr1",1,470535);
        GenomeLoc newLoc;
        while (i.hasNext()) {
            i.next();
            countOfReads++;
        }
        reader.close();
        while (last.getStart() < seq.getSequenceDictionary().getSequence("chr1").getSequenceLength()) {
            reader = new SAMFileReader(bam, true);
            long stop = (last.getStop() >= seqLength) ? seqLength : last.getStop() + 470535;
            newLoc = GenomeLocParser.createGenomeLoc(last.getContigIndex(),last.getStart()+470535,stop);
            reader.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
            i = reader.queryOverlapping(newLoc.getContig(),(int)newLoc.getStart(),(int)newLoc.getStop());
            IntervalOverlapIterator iter = new IntervalOverlapIterator(StingSAMIteratorAdapter.adapt(null, i),last,false);
            while(iter.hasNext()) {
                countOfReads++;
                iter.next();
            }
            last = newLoc;
            reader.close();

        }
        Assert.assertEquals(chromosomeOneReadCount,countOfReads);
    }
}
