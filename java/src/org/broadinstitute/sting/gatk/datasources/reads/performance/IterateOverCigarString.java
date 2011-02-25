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

package org.broadinstitute.sting.gatk.datasources.reads.performance;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Feb 25, 2011
 * Time: 10:16:53 AM
 * To change this template use File | Settings | File Templates.
 */
class IterateOverCigarString extends ReadProcessor {
    private long matchMismatches;
    private long insertions;
    private long deletions;
    private long others;

    public IterateOverCigarString(final BAMProcessingPerformanceMeter performanceMeter) {
        super(performanceMeter);
    }

    @Override
    public String getTestName() { return "iterator over cigar string"; }
    public void processRead(final SAMRecord read) {
        Cigar cigar = read.getCigar();
        for(CigarElement cigarElement: cigar.getCigarElements()) {
            int elementSize = cigarElement.getLength();
            while(elementSize > 0) {
                switch(cigarElement.getOperator()) {
                    case M: matchMismatches++; break;
                    case I: insertions++; break;
                    case D: deletions++; break;
                    default: others++; break;
                }
                elementSize--;
            }
        }
    }
}
