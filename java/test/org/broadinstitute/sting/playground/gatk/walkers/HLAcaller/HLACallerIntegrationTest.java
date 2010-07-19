/*
 * Copyright (c) 2010.
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

package org.broadinstitute.sting.playground.gatk.walkers.HLAcaller;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class HLACallerIntegrationTest extends WalkerTest {

    private static final String intervals = validationDataLocation + "HLA_EXONS.intervals";


    @Test
    public void testFindClosestHLA() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T FindClosestHLA -I " + validationDataLocation + "NA12878.HISEQ.HLA.bam -R " + oneKGLocation + "reference/human_b36_both.fasta -L " + intervals + " -useInterval " + intervals + " -HLAdictionary " + validationDataLocation + "HLA_DICTIONARY.txt -PolymorphicSites " + validationDataLocation + "HLA_POLYMORPHIC_SITES.txt -o %s", 1,
                Arrays.asList("a49b6f54a4585d1dd958c55a5523427d"));
        executeTest("test FindClosestHLA", spec);
    }

    @Test
    public void testCalculateBaseLikelihoods() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T CalculateBaseLikelihoods -I " + validationDataLocation + "NA12878.HISEQ.HLA.bam -R " + oneKGLocation + "reference/human_b36_both.fasta -L " + intervals + " -filter " + validationDataLocation + "HLA_HISEQ.filter -maxAllowedMismatches 6 -minRequiredMatches 0 -o %s", 1,
                Arrays.asList("98e64882f93bf7550457bee4182caab6"));
        executeTest("test CalculateBaseLikelihoods", spec);
    }

    @Test
    public void testHLACaller() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T HLACaller -noVerbose -I " + validationDataLocation + "NA12878.HISEQ.HLA.bam -R " + oneKGLocation + "reference/human_b36_both.fasta -L " + intervals + " -useInterval " + intervals + " -HLAdictionary " + validationDataLocation + "HLA_DICTIONARY.txt -filter " + validationDataLocation + "HLA_HISEQ.filter -maxAllowedMismatches 6 -minRequiredMatches 5 -HLAfrequencies " + validationDataLocation + "HLA_FREQUENCIES.txt -bl " + validationDataLocation + "HLA_HISEQ.baselikelihoods -o %s", 1,
                Arrays.asList("f9931b378bde213e71fca6ecaa24b48b"));
        executeTest("test HLACaller", spec);
    }
}