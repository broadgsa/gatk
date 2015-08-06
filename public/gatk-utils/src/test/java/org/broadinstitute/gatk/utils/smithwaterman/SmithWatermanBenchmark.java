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

package org.broadinstitute.gatk.utils.smithwaterman;

import com.google.caliper.Param;
import com.google.caliper.SimpleBenchmark;
import org.broadinstitute.gatk.utils.Utils;

/**
 * Caliper microbenchmark of parsing a VCF file
 */
public class SmithWatermanBenchmark extends SimpleBenchmark {

    @Param({"Original"})
    String version; // set automatically by framework

    @Param({"10", "50", "100", "500"})
    int sizeOfMiddleRegion; // set automatically by framework

    @Param({"10", "50", "100", "500"})
    int sizeOfEndRegions; // set automatically by framework

    String refString;
    String hapString;

    @Override protected void setUp() {
        final StringBuilder ref = new StringBuilder();
        final StringBuilder hap = new StringBuilder();

        ref.append(Utils.dupString('A', sizeOfEndRegions));
        hap.append(Utils.dupString('A', sizeOfEndRegions));

        // introduce a SNP
        ref.append("X");
        hap.append("Y");

        ref.append(Utils.dupString('A', sizeOfMiddleRegion));
        hap.append(Utils.dupString('A', sizeOfMiddleRegion));

        // introduce a SNP
        ref.append("X");
        hap.append("Y");

        ref.append(Utils.dupString('A', sizeOfEndRegions));
        hap.append(Utils.dupString('A', sizeOfEndRegions));

        refString = ref.toString();
        hapString = hap.toString();
    }

    public void timeSW(int rep) {
        for ( int i = 0; i < rep; i++ ) {
            final SmithWaterman sw;
            if ( version.equals("Greedy") )
                throw new IllegalArgumentException("Unsupported implementation");
            sw = new SWPairwiseAlignment(refString.getBytes(), hapString.getBytes());
            sw.getCigar();
        }
    }

    public static void main(String[] args) {
        com.google.caliper.Runner.main(SmithWatermanBenchmark.class, args);
    }
}
