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

package org.broadinstitute.gatk.engine.datasources.reads;

import com.google.caliper.Param;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;

import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Apr 22, 2011
 * Time: 4:01:23 PM
 * To change this template use File | Settings | File Templates.
 */
public class TheoreticalMinimaBenchmark extends ReadProcessingBenchmark {
    @Param
    private String bamFile;

    @Param
    private Integer maxReads;

    @Override
    public String getBAMFile() { return bamFile; }

    @Override
    public Integer getMaxReads() { return maxReads; }

    public void timeIterateOverEachBase(int reps) {
        System.out.printf("Processing " + inputFile);
        for(int i = 0; i < reps; i++) {
            SAMFileReader reader = new SAMFileReader(inputFile);
            CloseableIterator<SAMRecord> iterator = reader.iterator();

            long As=0,Cs=0,Gs=0,Ts=0;
            while(iterator.hasNext()) {
                SAMRecord read = iterator.next();
                for(byte base: read.getReadBases()) {
                    switch(base) {
                        case 'A': As++; break;
                        case 'C': Cs++; break;
                        case 'G': Gs++; break;
                        case 'T': Ts++; break;
                    }
                }
            }
            System.out.printf("As = %d; Cs = %d; Gs = %d; Ts = %d; total = %d%n",As,Cs,Gs,Ts,As+Cs+Gs+Ts);
            iterator.close();
            reader.close();
        }
    }

    public void timeIterateOverCigarString(int reps) {
        for(int i = 0; i < reps; i++) {
            long matchMismatches = 0;
            long insertions = 0;
            long deletions = 0;
            long others = 0;

            SAMFileReader reader = new SAMFileReader(inputFile);
            CloseableIterator<SAMRecord> iterator = reader.iterator();
            while(iterator.hasNext()) {
                SAMRecord read = iterator.next();

                Cigar cigar = read.getCigar();
                for(CigarElement cigarElement: cigar.getCigarElements()) {
                    int elementSize = cigarElement.getLength();
                    while(elementSize > 0) {
                        switch(cigarElement.getOperator()) {
                            case M: case EQ: case X: matchMismatches++; break;
                            case I: insertions++; break;
                            case D: deletions++; break;
                            default: others++; break;
                        }
                        elementSize--;
                    }
                }
            }
            System.out.printf("Ms = %d; Is = %d; Ds = %d; others = %d; total = %d%n",matchMismatches,insertions,deletions,others,matchMismatches+insertions+deletions+others);

            iterator.close();
            reader.close();
        }
    }

}
