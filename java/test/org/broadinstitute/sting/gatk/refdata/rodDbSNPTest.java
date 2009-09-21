package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.*;


/**
 * @author aaron
 *         <p/>
 *         Class rodDbSNPTest
 *         <p/>
 *         A descriptions should go here. Blame aaron if it's missing.
 */
public class rodDbSNPTest extends BaseTest {
    private static IndexedFastaSequenceFile seq;

    @BeforeClass
    public static void beforeTests() {
        try {
            seq = new IndexedFastaSequenceFile(new File("/broad/1KG/reference/human_b36_both.fasta"));
        } catch (FileNotFoundException e) {
            throw new StingException("unable to load the sequence dictionary");
        }
        GenomeLocParser.setupRefContigOrdering(seq);

    }

    public BufferedReader openFile(String filename) {
        try {
            return new BufferedReader(new FileReader(filename));
        } catch (FileNotFoundException e) {
            throw new StingException("Couldn't open file " + filename);
        }

    }

    @Test
    // count the number of SNP's between 10M and 11M on chr1 in dbSNP
    public void testDBSNPInput() {
        BufferedReader stream = openFile("/humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod");
        int snpCount = 0;
        int indelCount = 0;
        try {
            String line = stream.readLine();
            rodDbSNP rod = new rodDbSNP("db");
            boolean stop = false;
            while (line != null && !stop) {
                rod.parseLine(null,line.split("\t"));
                rodDbSNP var = (rodDbSNP)rod;
                if (rod.isSNP()) {
                    // quick check, if we're not triallelic, make sure the ref is right
                    if (var.getRefSnpFWD() == var.refBases.charAt(0) || var.getAltSnpFWD() == var.refBases.charAt(0))
                        // also make sure the ref is a single character
                        if (var.refBases.length() == 1)
                            Assert.assertTrue(var.refBases.charAt(0)==var.getRefSnpFWD());
                    if (var.getLocation().getContig().equals("1") &&
                            var.getLocation().getStart() >= 10000000 &&
                            var.getLocation().getStart() <= 11000000) {
                        if (var.isSNP()) {
                            snpCount++;
                        }
                        
                    }
                }
                if (rod.isIndel())
                    indelCount++;

                stop = (var.getLocation().getContig().equals("1") && var.getLocation().getStart() > 11000000);
                line = stream.readLine();
            }
            Assert.assertEquals(3615,snpCount);
            Assert.assertEquals(9902,indelCount);

        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }


}
