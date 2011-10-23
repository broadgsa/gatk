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

package org.broadinstitute.sting.utils.pileup;

import com.google.caliper.Param;
import com.google.caliper.SimpleBenchmark;
import com.google.caliper.runner.CaliperMain;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;

import java.util.*;

/**
 * Caliper microbenchmark of fragment pileup
 */
public class FragmentPileupBenchmark extends SimpleBenchmark {
    final int N_PILEUPS_TO_GENERATE = 100;
    List<ReadBackedPileup> pileups = new ArrayList<ReadBackedPileup>(N_PILEUPS_TO_GENERATE);

    @Param({"10", "100", "1000"}) // , "10000"})
    int pileupSize; // set automatically by framework

    @Param({"150", "400"})
    int insertSize; // set automatically by framework

    @Override protected void setUp() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        GenomeLocParser genomeLocParser;
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
        final int pos = 50;
        GenomeLoc loc = genomeLocParser.createGenomeLoc("chr1", pos);

        final Random ran = new Random();
        final int readLen = 100;
        final boolean leftIsFirst = true;
        final boolean leftIsNegative = false;
        final int insertSizeVariation = insertSize / 10;

        for ( int pileupN = 0; pileupN < N_PILEUPS_TO_GENERATE; pileupN++ ) {
            List<PileupElement> pileupElements = new ArrayList<PileupElement>();
            for ( int i = 0; i < pileupSize / 2; i++ ) {
                final String readName = "read" + i;
                final int leftStart = new Random().nextInt(49) + 1;
                final int fragmentSize = (int)(ran.nextGaussian() * insertSizeVariation + insertSize);
                final int rightStart = leftStart + fragmentSize - readLen;

                if ( rightStart <= 0 ) continue;

                List<SAMRecord> pair = FragmentPileupUnitTest.createPair(header, readName, readLen, leftStart, rightStart, leftIsFirst, leftIsNegative);
                SAMRecord left = pair.get(0);
                SAMRecord right = pair.get(1);

                pileupElements.add(new PileupElement(left, pos - leftStart));

                if ( pos >= right.getAlignmentStart() && pos <= right.getAlignmentEnd() ) {
                    pileupElements.add(new PileupElement(right, pos - rightStart));
                }
            }

            Collections.sort(pileupElements);
            pileups.add(new ReadBackedPileupImpl(loc, pileupElements));
        }
    }

    public void timeNaiveNameMatch(int rep) {
        int nFrags = 0;
        for ( int i = 0; i < rep; i++ ) {
            for ( ReadBackedPileup rbp : pileups )
                nFrags += new FragmentPileup(rbp, true).getTwoReadPileup().size();
        }
    }

    public void timeFastNameMatch(int rep) {
        int nFrags = 0;
        for ( int i = 0; i < rep; i++ )
            for ( ReadBackedPileup rbp : pileups )
                nFrags += new FragmentPileup(rbp, false).getTwoReadPileup().size();
    }

    public static void main(String[] args) {
        CaliperMain.main(FragmentPileupBenchmark.class, args);
    }
}
