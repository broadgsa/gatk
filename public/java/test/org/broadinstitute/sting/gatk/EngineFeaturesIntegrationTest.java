/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 *
 */
public class EngineFeaturesIntegrationTest extends WalkerTest {
    private void testBadRODBindingInput(String type, String name, Class c) {
        WalkerTestSpec spec = new WalkerTestSpec("-T SelectVariants -L 1:1 --variant:variant," + type + " "
                + b37dbSNP132 + " -R " + b37KGReference + " -o %s",
                1, c);
        executeTest(name, spec);
    }

    @Test() private void testBadRODBindingInputType1() {
        testBadRODBindingInput("beagle", "BEAGLE input to VCF expecting walker", UserException.BadArgumentValue.class);
    }

    @Test() private void testBadRODBindingInputType3() {
        testBadRODBindingInput("bed", "Bed input to VCF expecting walker", UserException.BadArgumentValue.class);
    }

    @Test() private void testBadRODBindingInputTypeUnknownType() {
        testBadRODBindingInput("bedXXX", "Unknown input to VCF expecting walker", UserException.UnknownTribbleType.class);
    }

    private void testMissingFile(String name, String missingBinding) {
        WalkerTestSpec spec = new WalkerTestSpec(missingBinding + " -R " + b37KGReference + " -o %s",
                1, UserException.CouldNotReadInputFile.class);
        executeTest(name, spec);
    }

    @Test() private void testMissingBAMnt1() {
        testMissingFile("missing BAM", "-T UnifiedGenotyper -I missing.bam -nt 1");
    }
    @Test() private void testMissingBAMnt4() {
        testMissingFile("missing BAM", "-T UnifiedGenotyper -I missing.bam -nt 4");
    }
    @Test() private void testMissingVCF() {
        testMissingFile("missing VCF", "-T SelectVariants -V missing.vcf");
    }
    @Test() private void testMissingInterval() {
        testMissingFile("missing interval", "-T UnifiedGenotyper -L missing.interval_list -I " + b37GoodBAM);
    }


    // --------------------------------------------------------------------------------
    //
    // Test that our exceptions are coming back as we expect
    //
    // --------------------------------------------------------------------------------

    private class EngineErrorHandlingTestProvider extends TestDataProvider {
        final Class expectedException;
        final boolean multiThreaded;
        final int iterationsToTest;

        public EngineErrorHandlingTestProvider(Class exceptedException, final boolean multiThreaded) {
            super(EngineErrorHandlingTestProvider.class);
            this.expectedException = exceptedException;
            this.multiThreaded = multiThreaded;
            this.iterationsToTest = multiThreaded ? 10 : 1;
            setName(String.format("Engine error handling: expected %s, is-multithreaded %b", exceptedException, multiThreaded));
        }
    }

    @DataProvider(name = "EngineErrorHandlingTestProvider")
    public Object[][] makeEngineErrorHandlingTestProvider() {
        for ( final boolean multiThreaded : Arrays.asList(true, false)) {
            new EngineErrorHandlingTestProvider(NullPointerException.class, multiThreaded);
            new EngineErrorHandlingTestProvider(UserException.class, multiThreaded);
            new EngineErrorHandlingTestProvider(ReviewedStingException.class, multiThreaded);
        }

        return EngineErrorHandlingTestProvider.getTests(EngineErrorHandlingTestProvider.class);
    }

    //
    // Loop over errors to throw, make sure they are the errors we get back from the engine, regardless of NT type
    //
    @Test(dataProvider = "EngineErrorHandlingTestProvider")
    public void testEngineErrorHandlingTestProvider(EngineErrorHandlingTestProvider cfg) {
        for ( int i = 0; i < cfg.iterationsToTest; i++ ) {
            final String root = "-T ErrorThrowing -R " + b37KGReference;
            final String args = root + (cfg.multiThreaded ? " -nt 2" : "") + " -E " + cfg.expectedException.getSimpleName();
            WalkerTestSpec spec = new WalkerTestSpec(args, 0, cfg.expectedException);
            executeTest(cfg.toString(), spec);
        }
    }
}