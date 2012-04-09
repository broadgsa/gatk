package org.broadinstitute.sting.utils.recalibration;

import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.utils.NGSPlatform;
import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Unit tests for on-the-fly recalibration.
 *
 * @author Mauricio Carneiro
 * @since 3/16/12
 */
public class BaseRecalibrationUnitTest {

    @Test(enabled=false)
    public void testReadingReport() {
        File csv = new File("public/testdata/exampleGATKREPORT.grp");
        BaseRecalibration baseRecalibration = new BaseRecalibration(csv, -1);
        GATKSAMRecord read = ReadUtils.createRandomRead(1000);
        read.setReadGroup(new GATKSAMReadGroupRecord(new SAMReadGroupRecord("exampleBAM.bam.bam"), NGSPlatform.ILLUMINA));
        baseRecalibration.recalibrateRead(read);
        System.out.println("Success");
    }
}
