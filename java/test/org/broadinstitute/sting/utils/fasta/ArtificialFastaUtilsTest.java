package org.broadinstitute.sting.utils.fasta;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;


/**
 * @author aaron
 *         <p/>
 *         Class ArtificialFastaUtilsTest
 *         <p/>
 *         test out the ArtificialFastaUtils functionality
 */
public class ArtificialFastaUtilsTest extends BaseTest {

    /** generate a fake fasta */
    @Test
    public void testFastaGeneration() {
        List<String> names = new ArrayList<String>();
        List<Integer> sizes = new ArrayList<Integer>();

        for (int x = 0; x < 5; x++) {
            sizes.add(1000);
            names.add("chr" + (x+1));
        }
        File temp = new File("tempFileFasta.fasta");
        ArtificialFastaUtils.createArtificialFasta(temp.getName(),names,sizes,ArtificialFastaUtils.BASE_PATTERN.ALL_A);

        // using the fasta sequence file to test, in reality we should use the indexed version
        FastaSequenceFile2 fasta = new FastaSequenceFile2(temp);

        Assert.assertEquals(5,fasta.getSequenceDictionary().getSequences().size());

        ArtificialSAMUtils.createArtificialBamFile("tempFileBAM.bam",5,1,1000,600);
        //temp.delete();
    }
}
