/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.gatk.tools;

import org.apache.commons.lang.StringUtils;
import htsjdk.tribble.AbstractFeatureReader;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.MD5DB;
import org.broadinstitute.gatk.utils.MD5Mismatch;
import org.broadinstitute.gatk.utils.runtime.ProcessController;
import org.broadinstitute.gatk.utils.runtime.ProcessSettings;
import org.broadinstitute.gatk.utils.runtime.RuntimeUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class CatVariantsIntegrationTest {
    private final MD5DB md5db = new MD5DB();
    private final File CatVariantsDir = new File(BaseTest.privateTestDir, "CatVariants");
    private final File CatVariantsVcf1 = new File(CatVariantsDir, "CatVariantsTest1.vcf");
    private final File CatVariantsVcf2 = new File(CatVariantsDir, "CatVariantsTest2.vcf");
    private final File CatVariantsBcf1 = new File(CatVariantsDir, "CatVariantsTest1.bcf");
    private final File CatVariantsBcf2 = new File(CatVariantsDir, "CatVariantsTest2.bcf");

    private class CatVariantsTestProvider extends BaseTest.TestDataProvider {
        private final File file1;
        private final File file2;
        public final File outputFile;
        public final String md5;

        private CatVariantsTestProvider(final File file1, final File file2, final File outputFile, final String md5) {
            super(CatVariantsTestProvider.class);

            this.file1 = file1;
            this.file2 = file2;
            this.outputFile = outputFile;
            this.md5 = md5;
        }

        public final String getCmdLine() {
            return String.format("java -cp %s %s -R %s -V %s -V %s -out %s",
                    StringUtils.join(RuntimeUtils.getAbsoluteClassPaths(), File.pathSeparatorChar),
                    CatVariants.class.getCanonicalName(), BaseTest.b37KGReference, file1, file2, outputFile);
        }

        public String toString() {
            return "CatVariantsTestProvider " + outputFile;
        }
    }

    @DataProvider(name = "ExtensionsTest")
    public Object[][] makeExtensionsTestProvider() {
        final File catVariantsTempList1 = BaseTest.createTempListFile("CatVariantsTest1", CatVariantsVcf1.getAbsolutePath());
        final File catVariantsTempList2 = BaseTest.createTempListFile("CatVariantsTest2", CatVariantsVcf2.getAbsolutePath());

        new CatVariantsTestProvider(CatVariantsVcf1, CatVariantsVcf2, BaseTest.createTempFile("CatVariantsTest", ".vcf"), "d0d81eb7fd3905256c4ac7c0fc480094");
        new CatVariantsTestProvider(CatVariantsBcf1, CatVariantsBcf2, BaseTest.createTempFile("CatVariantsTest", ".bcf"), "6a57fcbbf3cae490896d13a288670d83");

        for (String extension : AbstractFeatureReader.BLOCK_COMPRESSED_EXTENSIONS) {
            final File file1 = new File(CatVariantsDir, "CatVariantsTest1.vcf" + extension);
            final File file2 = new File(CatVariantsDir, "CatVariantsTest2.vcf" + extension);
            final File outputFile = BaseTest.createTempFile("CatVariantsTest", ".vcf" + extension);
            new CatVariantsTestProvider(file1, file2, outputFile, "33f728ac5c70ce2994f3619a27f47088");
        }

        //Test list parsing functionality
        new CatVariantsTestProvider(catVariantsTempList1, CatVariantsVcf2, BaseTest.createTempFile("CatVariantsTest", ".vcf"), "d0d81eb7fd3905256c4ac7c0fc480094");
        new CatVariantsTestProvider(CatVariantsVcf1, catVariantsTempList2, BaseTest.createTempFile("CatVariantsTest", ".vcf"), "d0d81eb7fd3905256c4ac7c0fc480094");
        new CatVariantsTestProvider(catVariantsTempList1, catVariantsTempList2, BaseTest.createTempFile("CatVariantsTest", ".vcf"), "d0d81eb7fd3905256c4ac7c0fc480094");

        return CatVariantsTestProvider.getTests(CatVariantsTestProvider.class);
    }

    @Test(dataProvider = "ExtensionsTest")
    public void testExtensions(final CatVariantsTestProvider cfg) throws IOException {

        ProcessController pc = ProcessController.getThreadLocal();
        ProcessSettings ps = new ProcessSettings(cfg.getCmdLine().split("\\s+"));
        pc.execAndCheck(ps);

        MD5DB.MD5Match result = md5db.testFileMD5("testExtensions", "CatVariantsTestProvider", cfg.outputFile, cfg.md5, false);
        if(result.failed) {
            final MD5Mismatch failure = new MD5Mismatch(result.actualMD5, result.expectedMD5, result.diffEngineOutput);
            Assert.fail(failure.toString());
        }
    }

    @Test(expectedExceptions = IOException.class)
    public void testMismatchedExtensions1() throws IOException {

        String cmdLine = String.format("java -cp %s %s -R %s -V %s -V %s -out %s",
                StringUtils.join(RuntimeUtils.getAbsoluteClassPaths(), File.pathSeparatorChar),
                CatVariants.class.getCanonicalName(),
                BaseTest.b37KGReference,
                CatVariantsVcf1,
                CatVariantsVcf2,
                BaseTest.createTempFile("CatVariantsTest", ".bcf"));

        ProcessController pc = ProcessController.getThreadLocal();
        ProcessSettings ps = new ProcessSettings(cmdLine.split("\\s+"));
        pc.execAndCheck(ps);
    }

    @Test(expectedExceptions = IOException.class)
    public void testMismatchedExtensions2() throws IOException {

        String cmdLine = String.format("java -cp %s %s -R %s -V %s -V %s -out %s",
                StringUtils.join(RuntimeUtils.getAbsoluteClassPaths(), File.pathSeparatorChar),
                CatVariants.class.getCanonicalName(),
                BaseTest.b37KGReference,
                CatVariantsVcf1,
                CatVariantsBcf2,
                BaseTest.createTempFile("CatVariantsTest", ".vcf"));

        ProcessController pc = ProcessController.getThreadLocal();
        ProcessSettings ps = new ProcessSettings(cmdLine.split("\\s+"));
        pc.execAndCheck(ps);
    }
}