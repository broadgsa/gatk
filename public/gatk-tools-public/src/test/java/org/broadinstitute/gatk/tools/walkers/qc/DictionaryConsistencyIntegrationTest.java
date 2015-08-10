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

package org.broadinstitute.gatk.tools.walkers.qc;

import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.annotations.Test;

/**
 *
 */
public class DictionaryConsistencyIntegrationTest extends WalkerTest {
    private static final String callsHG18 = BaseTest.validationDataLocation + "HiSeq.WGS.cleaned.ug.snpfiltered.indelfiltered.optimized.cut.vcf";
    private static final String callsB36  = BaseTest.validationDataLocation + "lowpass.N3.chr1.raw.vcf";
    private static final String lexHG18 = BaseTest.validationDataLocation + "lexFasta/lex.hg18.fasta";
    private static final String b36BAM = BaseTest.validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam";
    private static final String hg18BAM = BaseTest.validationDataLocation + "NA12878.WEx.downsampled20x.bam";

    // testing hg18 calls
    @Test public void fail1() { executeTest("b36xhg18", testVCF(b36KGReference, callsHG18)); }
    @Test public void fail2() { executeTest("b37xhg18", testVCF(b37KGReference, callsHG18)); }
    @Test public void fail3() { executeTest("hg19xhg18", testVCF(hg19Reference, callsHG18)); }
    @Test public void fail4() { executeTest("hg18lex-v-hg18", testVCF(lexHG18, callsHG18, UserException.LexicographicallySortedSequenceDictionary.class)); }

    // testing b36 calls
    @Test public void fail5() { executeTest("b36xb36", testVCF(hg18Reference, callsB36)); }
    @Test public void fail6() { executeTest("b37xb36", testVCF(b37KGReference, callsB36)); }
    @Test public void fail7() { executeTest("hg19xb36", testVCF(hg19Reference, callsB36)); }
    @Test public void fail8() { executeTest("hg18lex-v-b36", testVCF(lexHG18, callsB36, UserException.LexicographicallySortedSequenceDictionary.class)); }

    private WalkerTest.WalkerTestSpec testVCF(String ref, String vcf) {
        return testVCF(ref, vcf, UserException.IncompatibleSequenceDictionaries.class);
    }

    private WalkerTest.WalkerTestSpec testVCF(String ref, String vcf, Class c) {
        return new WalkerTest.WalkerTestSpec("-T VariantsToTable -M 10 --variant:vcf "
                + vcf + " -F POS,CHROM -R "
                + ref +  " -o %s",
                1, c);

    }

    @Test public void failBAM1() { executeTest("b36bam-v-b37",     testBAM(b37KGReference, b36BAM, "1",    UserException.IncompatibleSequenceDictionaries.class)); }
    @Test public void failBAM2() { executeTest("b36bam-v-hg18",    testBAM(hg18Reference,  b36BAM, "chr1", UserException.IncompatibleSequenceDictionaries.class)); }
    @Test public void failBAM3() { executeTest("b36bam-v-hg19",    testBAM(hg19Reference,  b36BAM, "1",    UserException.IncompatibleSequenceDictionaries.class)); }
    @Test public void failBAM4() { executeTest("b36bam-v-lexhg18", testBAM(lexHG18,        b36BAM, "chr1", UserException.LexicographicallySortedSequenceDictionary.class)); }

    @Test public void failBAM5() { executeTest("hg18bam-v-b36",     testBAM(b36KGReference, hg18BAM, "1",    UserException.IncompatibleSequenceDictionaries.class)); }
    @Test public void failBAM6() { executeTest("hg18bam-v-b37",     testBAM(b37KGReference, hg18BAM, "1",    UserException.IncompatibleSequenceDictionaries.class)); }
    @Test public void failBAM7() { executeTest("hg18bam-v-hg19",    testBAM(hg19Reference,  hg18BAM, "1",    UserException.IncompatibleSequenceDictionaries.class)); }
    @Test public void failBAM8() { executeTest("hg18bam-v-lexhg18", testBAM(lexHG18,        hg18BAM, "chr1", UserException.LexicographicallySortedSequenceDictionary.class)); }

    private WalkerTest.WalkerTestSpec testBAM(String ref, String bam, String contig, Class c) {
        return new WalkerTest.WalkerTestSpec("-T PrintReads -I " + bam + " -R " + ref +  " -L " + contig + ":10,000,000-11,000,000 -o %s",
                1, c);

    }
}
