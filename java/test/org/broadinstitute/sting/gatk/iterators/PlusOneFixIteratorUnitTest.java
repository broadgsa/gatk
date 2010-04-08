package org.broadinstitute.sting.gatk.iterators;

import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.BaseTest;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.Assert;
import static org.junit.Assert.fail;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Iterator;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: aaronmckenna
 * Date: Nov 4, 2009
 * Time: 11:02:24 PM
 */
public class PlusOneFixIteratorUnitTest extends BaseTest {
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

    /**
     * there's a read starting at chr1:1108664, make our interval go from chr1 to 1108664, apply the iterator,
     * and we shouldn't see the read.
     */
    @Test
    public void testReadAtEndOfInterval() {
        int countOfReads = 0;
        SAMFileReader reader = new SAMFileReader(bam, true);
        final int size = 1108663;
        reader.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);

        GenomeLoc last = GenomeLocParser.createGenomeLoc("chr1", 1, size);
        Iterator<SAMRecord> i = new PlusOneFixIterator(last, StingSAMIteratorAdapter.adapt(null, reader.queryOverlapping("chr1", 1, 1108663)));
        //Iterator<SAMRecord> i = reader.queryOverlapping("chr1", 1, size);
        while (i.hasNext()) {
            SAMRecord rec = i.next();
            countOfReads++;
        }
        Assert.assertEquals(2, countOfReads);
    }

    /**
     * there are 4296 reads in chr1, make sure we see them all.
     *
     * chr1 length = 247249719 (and no reads start at the last position
     */
    @Test
    public void testAllReads() {
        int countOfReads = 0;
        SAMFileReader reader = new SAMFileReader(bam, true);
        reader.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
        final int size = 247249719;
        GenomeLoc last = GenomeLocParser.createGenomeLoc("chr1", 1, size);
        Iterator<SAMRecord> i = new PlusOneFixIterator(last, StingSAMIteratorAdapter.adapt(null, reader.queryOverlapping("chr1", 1, size)));
        //Iterator<SAMRecord> i = reader.queryOverlapping("chr1", 1, size);
        while (i.hasNext()) {
            SAMRecord rec = i.next();
            countOfReads++;
        }
        Assert.assertEquals(885, countOfReads);
    }

    /**
     * there are 4296 reads in chr1, make sure we see them all.
     *
     * chr1 length = 247249719 (and no reads start at the last position
     */
    @Test
    public void testAllReadsSharded() {
        int countOfReads = 0;

        final int size = 247249719;
        int currentStop = 0;
        int incr = 100000;
        GenomeLoc last = GenomeLocParser.createGenomeLoc("chr1", 1, size);
        while (currentStop < size) {
            SAMFileReader reader = new SAMFileReader(bam, true);
        reader.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
            int lastStop = (currentStop > 1) ? currentStop: 1;
            if (currentStop + incr < size)
                currentStop += incr;
            else
                currentStop = size;
            Iterator<SAMRecord> i = new PlusOneFixIterator(last, StingSAMIteratorAdapter.adapt(null, reader.queryOverlapping("chr1", lastStop, currentStop)));

            while (i.hasNext()) {
                SAMRecord rec = i.next();
                countOfReads++;
            }            
            reader.close();
        }
        Assert.assertEquals(885, countOfReads);
    }
}

