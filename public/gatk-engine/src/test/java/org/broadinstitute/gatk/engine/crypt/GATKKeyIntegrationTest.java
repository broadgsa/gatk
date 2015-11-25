/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.crypt;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.broadinstitute.gatk.engine.phonehome.GATKRunReport;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

public class GATKKeyIntegrationTest extends WalkerTest {

    public static final String BASE_COMMAND = String.format("-T TestPrintReadsWalker -R %s -I %s -o %%s",
                                                            publicTestDir + "exampleFASTA.fasta",
                                                            publicTestDir + "exampleBAM.bam");
    public static final String MD5_UPON_SUCCESSFUL_RUN = "462656ec9632f8c21ee534d35093c3f8";


    private void runGATKKeyTest ( String testName, String etArg, String keyArg, Class expectedException, String md5 ) {
        String command = BASE_COMMAND + String.format(" %s %s", etArg, keyArg);

        WalkerTestSpec spec = expectedException != null ?
                              new WalkerTestSpec(command, 1, expectedException) :
                              new WalkerTestSpec(command, 1, Arrays.asList(md5));

        spec.disableImplicitArgs(); // Turn off automatic inclusion of -et/-K args by WalkerTest
        executeTest(testName, spec);
    }

    @Test
    public void testValidKeyNoET() {
        runGATKKeyTest("testValidKeyNoET",
                       "-et " + GATKRunReport.PhoneHomeOption.NO_ET,
                       "-K " + keysDataLocation + "valid.key",
                       null,
                       MD5_UPON_SUCCESSFUL_RUN);
    }

    @Test
    public void testValidKeyETStdout() {
        runGATKKeyTest("testValidKeyETStdout",
                       "-et " + GATKRunReport.PhoneHomeOption.STDOUT,
                       "-K " + keysDataLocation + "valid.key",
                       null,
                       MD5_UPON_SUCCESSFUL_RUN);
    }

    @Test
    public void testValidKeyETStandard() {
        runGATKKeyTest("testValidKeyETStandard",
                       "",
                       "-K " + keysDataLocation + "valid.key",
                       null,
                       MD5_UPON_SUCCESSFUL_RUN);
    }

    @Test
    public void testNoKeyNoET() {
        runGATKKeyTest("testNoKeyNoET",
                       "-et " + GATKRunReport.PhoneHomeOption.NO_ET,
                       "",
                       UserException.class,
                       null);
    }

    @Test
    public void testNoKeyETStdout() {
        runGATKKeyTest("testNoKeyETStdout",
                       "-et " + GATKRunReport.PhoneHomeOption.STDOUT,
                       "",
                       UserException.class,
                       null);
    }

    @Test
    public void testNoKeyETStandard() {
        runGATKKeyTest("testNoKeyETStandard",
                       "",
                       "",
                       null,
                       MD5_UPON_SUCCESSFUL_RUN);
    }

    @Test
    public void testRevokedKey() {
        runGATKKeyTest("testRevokedKey",
                       "-et " + GATKRunReport.PhoneHomeOption.NO_ET,
                       "-K " + keysDataLocation + "revoked.key",
                       UserException.KeySignatureVerificationException.class,
                       null);
    }

    @DataProvider(name = "CorruptKeyTestData")
    public Object[][] corruptKeyDataProvider() {
        return new Object[][] {
            { "corrupt_empty.key",                  UserException.UnreadableKeyException.class },
            { "corrupt_single_byte_file.key",       UserException.UnreadableKeyException.class },
            { "corrupt_random_contents.key",        UserException.UnreadableKeyException.class },
            { "corrupt_single_byte_deletion.key",   UserException.UnreadableKeyException.class },
            { "corrupt_single_byte_insertion.key",  UserException.UnreadableKeyException.class },
            { "corrupt_single_byte_change.key",     UserException.UnreadableKeyException.class },
            { "corrupt_multi_byte_deletion.key",    UserException.UnreadableKeyException.class },
            { "corrupt_multi_byte_insertion.key",   UserException.UnreadableKeyException.class },
            { "corrupt_multi_byte_change.key",      UserException.UnreadableKeyException.class },
            { "corrupt_bad_isize_field.key",        UserException.UnreadableKeyException.class },
            { "corrupt_bad_crc.key",                UserException.UnreadableKeyException.class },
            { "corrupt_no_email_address.key",       UserException.UnreadableKeyException.class },
            { "corrupt_no_sectional_delimiter.key", UserException.UnreadableKeyException.class },
            { "corrupt_no_signature.key",           UserException.UnreadableKeyException.class },
            { "corrupt_bad_signature.key",          UserException.KeySignatureVerificationException.class },
            { "corrupt_non_gzipped_valid_key.key",  UserException.UnreadableKeyException.class }
        };
    }

    @Test(dataProvider = "CorruptKeyTestData")
    public void testCorruptKey ( String corruptKeyName, Class expectedException ) {
        runGATKKeyTest(String.format("testCorruptKey (%s)", corruptKeyName),
                       "-et " + GATKRunReport.PhoneHomeOption.NO_ET,
                       "-K " + keysDataLocation + corruptKeyName,
                       expectedException,
                       null);
    }

    @Test
    public void testCorruptButNonRequiredKey() {
        runGATKKeyTest("testCorruptButNonRequiredKey",
                       "",
                       "-K " + keysDataLocation + "corrupt_random_contents.key",
                       null,
                       MD5_UPON_SUCCESSFUL_RUN);
    }
}
