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

package org.broadinstitute.sting.tools;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.MD5DB;
import org.broadinstitute.sting.MD5Mismatch;
import org.broadinstitute.sting.utils.runtime.ProcessController;
import org.broadinstitute.sting.utils.runtime.ProcessSettings;
import org.broadinstitute.sting.utils.runtime.RuntimeUtils;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class CatVariantsIntegrationTest {
    private final MD5DB md5db = new MD5DB();
    private final File CatVariantsDir = new File(BaseTest.privateTestDir, "CatVariants");

    private class CatVariantsTestProvider extends BaseTest.TestDataProvider {
        private final File file1;
        private final File file2;
        public final File outputFile;
        public final String md5;

        private CatVariantsTestProvider(final String file1, final String file2, final File outputFile, final String md5) {
            super(CatVariantsTestProvider.class);

            this.file1 = new File(CatVariantsDir, file1);
            this.file2 = new File(CatVariantsDir, file2);
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
        new CatVariantsTestProvider("CatVariantsTest1.vcf", "CatVariantsTest2.vcf", BaseTest.createTempFile("CatVariantsTest", ".vcf"), "d0d81eb7fd3905256c4ac7c0fc480094");
        new CatVariantsTestProvider("CatVariantsTest1.bcf", "CatVariantsTest2.bcf", BaseTest.createTempFile("CatVariantsTest", ".bcf"), "6a57fcbbf3cae490896d13a288670d83");

        for (String extension : VariantContextWriterFactory.BLOCK_COMPRESSED_EXTENSIONS)
            new CatVariantsTestProvider("CatVariantsTest1.vcf" + extension, "CatVariantsTest2.vcf" + extension, BaseTest.createTempFile("CatVariantsTest", ".vcf" + extension), "33f728ac5c70ce2994f3619a27f47088");

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
                new File(CatVariantsDir, "CatVariantsTest1.vcf"),
                new File(CatVariantsDir, "CatVariantsTest2.vcf"),
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
                new File(CatVariantsDir, "CatVariantsTest1.vcf"),
                new File(CatVariantsDir, "CatVariantsTest2.bcf"),
                BaseTest.createTempFile("CatVariantsTest", ".vcf"));

        ProcessController pc = ProcessController.getThreadLocal();
        ProcessSettings ps = new ProcessSettings(cmdLine.split("\\s+"));
        pc.execAndCheck(ps);
    }
}