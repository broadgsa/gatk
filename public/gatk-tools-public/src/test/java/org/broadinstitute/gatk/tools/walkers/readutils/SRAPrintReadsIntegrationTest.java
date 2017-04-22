/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.readutils;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.testng.SkipException;
import org.testng.annotations.BeforeGroups;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;

import htsjdk.samtools.sra.SRAAccession;

import java.lang.reflect.Method;

public class SRAPrintReadsIntegrationTest extends WalkerTest {

    private static class PRTest {
        final String accession;
        final String args;
        final String md5;

        private PRTest(String accession, String args, String md5) {
            this.accession = accession;
            this.args = args;
            this.md5 = md5;
        }

        @Override
        public String toString() {
            return String.format("PRTest(accession='%s', args='%s')", accession, args);
        }
    }

    private static boolean canResolveNetworkAccession = false;
    private static String checkAccession = "SRR000123";
    private static final String REMOTE_ACCESSION_PATTERN = "^[SED]RR[0-9]{6,9}$";

    /**
     * Are the SRA native libraries loaded and initialized? Does the test accession have a valid name?
     */
    @BeforeGroups(groups = {"sra"})
    public final void checkIfCanResolve() {
        // Did SRA successfully load the native libraries and are fully initialized?
        if (!SRAAccession.isSupported()) {
            return;
        }

        // Is this is a valid SRA accession?
        canResolveNetworkAccession = SRAAccession.isValid(checkAccession);
    }

    /**
     * Are the SRA native libraries loaded and initialized?
     *
     * @throws SkipException if the SRA native libraries are loaded and initialized
     */
    @BeforeMethod
    public final void assertSRAIsSupported(){
        if(!SRAAccession.isSupported()){
            throw new SkipException("Skipping SRA Test because SRA native code is unavailable.");
        }
    }

    /**
     * Skip network SRA Test because cannot resolve remote SRA accession
     *
     * @param method    Provides information about, and access to, a single method on a class or interface
     * @param params    Method parameters
     * @throws SkipException if cannot resold an SRA accession
     */
    @BeforeMethod
    public  void skipIfCantResolve(Method method, Object[] params) {
        String accession = null;

        if (params.length > 0) {
            Object firstParam = params[0];
            if (firstParam instanceof PRTest) {
                accession = ((PRTest)firstParam).accession;
            }
        }

        if (accession != null &&
                accession.matches(REMOTE_ACCESSION_PATTERN) && !canResolveNetworkAccession) {
            throw new SkipException("Skipping network SRA Test because cannot resolve remote SRA accession '" +
                    checkAccession + "'.");
        }
    }

    @DataProvider(name = "PRTest")
    public Object[][] createPrintReadsTestData() {
        return new Object[][]{
                // Test with local SRA accessions
                {new PRTest(privateTestDir + "ERR1214757.sra", "",  "173ed87acc794a704aa000c6ab5d63a8")},
                {new PRTest(privateTestDir + "ERR1214757.sra", " -L NC_000001.10:1-50000",  "6bc055f028c49bcbca990857e57a6e4b")},
                {new PRTest(privateTestDir + "ERR1214757.sra", " -L NC_000001.10:500000-1000000",  "ab545064b2314971cfae7486ff74d779")},
                // Tests with remote SRA accessions
                {new PRTest("ERR1214757", "",  "173ed87acc794a704aa000c6ab5d63a8")},
                {new PRTest("ERR1214757", " -L NC_000001.10:1-50000",  "6bc055f028c49bcbca990857e57a6e4b")},
                {new PRTest("ERR1214757", " -L NC_000001.10:500000-1000000",  "ab545064b2314971cfae7486ff74d779")},
        };
    }

    @Test(groups = "sra", dataProvider = "PRTest")
    public void testPrintReads(PRTest params) {

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T PrintReads" +
                        " -R " + params.accession +
                        " -I " + params.accession +
                        params.args +
                        " --no_pg_tag" +
                        " -o %s",
                Collections.singletonList(params.md5));
        executeTest("testPrintReads-"+params.args, spec).getFirst();
    }

}
