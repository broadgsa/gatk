/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.alignment.bwa.BWAAligner;
import org.broadinstitute.sting.alignment.bwa.BWAConfiguration;
import org.broadinstitute.sting.alignment.bwa.c.BWACAligner;
import org.broadinstitute.sting.alignment.Alignment;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.utils.GATKException;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.StringUtil;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.SortedMap;

/**
 * A walker to experiment with fishing for reads in the GATK.  Has very limited utility in its current state.
 *
 * @author mhanna
 * @version 0.1
 */
public class TestReadFishingWalker extends ReadWalker<Integer,Long> {
    /**
     * An aligner for the small custom reference.
     */
    private BWAAligner aligner;

    @Output
    private PrintStream out;

    @Argument(fullName="indel_calls",shortName="ic",doc="Indel calls to use to derive custom references",required=true)
    private File indelCalls;

    @Argument(fullName="buffer_width",shortName="bw",doc="How much reference to extract around the given event",required=false)
    private int bufferWidth = 36;

    private SortedMap<GenomeLoc,BWAAligner> aligners = new TreeMap<GenomeLoc,BWAAligner>();

    @Override
    public void initialize() {
        long startTime = System.currentTimeMillis();
        int numAlignersCreated = 0;

        IndexedFastaSequenceFile referenceReader;
        FileInputStream indelCallInputStream;
        try {
            referenceReader = new IndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
            indelCallInputStream = new FileInputStream(indelCalls);
        }
        catch(IOException ex) {
            throw new UserException.CouldNotReadInputFile(indelCalls, ex);
        }

        Scanner indelCallReader = new Scanner(indelCallInputStream);

        while(indelCallReader.hasNext()) {
            String contig = indelCallReader.next();
            int eventPos = indelCallReader.nextInt();
            int eventLength = indelCallReader.nextInt();
            char type = indelCallReader.next().toUpperCase().charAt(0);
            byte[] bases = StringUtil.stringToBytes(indelCallReader.next());
            String sample = indelCallReader.next();

            byte[] revisedReference;
            int start,stop;

            if(type == 'D') {
                start = eventPos-eventLength-bufferWidth;
                stop = eventPos+eventLength+bufferWidth;
                int eventStart = eventPos - start + 1;
                int eventStop = eventStart + eventLength - 1;

                ReferenceSequence referenceSequence = referenceReader.getSubsequenceAt(contig,start,stop);
                revisedReference = new byte[(stop-start+1) - eventLength];

                System.arraycopy(referenceSequence.getBases(),0,revisedReference,0,eventStart);
                System.arraycopy(referenceSequence.getBases(),eventStop+1,revisedReference,eventStart,stop-start-eventStop);

            }
            else if(type == 'I') {
                start = eventPos-bufferWidth;
                stop = eventPos+bufferWidth;
                int eventStart = eventPos - start + 1;

                ReferenceSequence referenceSequence = referenceReader.getSubsequenceAt(contig,start,stop);
                revisedReference = new byte[(stop-start+1) + eventLength];

                System.arraycopy(referenceSequence.getBases(),0,revisedReference,0,bufferWidth+1);
                System.arraycopy(bases,0,revisedReference,eventStart,eventLength);
                System.arraycopy(referenceSequence.getBases(),eventStart,revisedReference,eventStart+eventLength,bufferWidth);
            }
            else
                throw new GATKException("Invalid indel type: " + type);

            aligners.put(GenomeLocParser.createGenomeLoc(contig,start,stop),new BWACAligner(revisedReference,new BWAConfiguration()));
            if(++numAlignersCreated % 100 == 0)
                out.printf("Created %d aligners in %dms%n",++numAlignersCreated,System.currentTimeMillis()-startTime);
        }
    }

    @Override
    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        Alignment bestAlignment = aligner.getBestAlignment(read.getReadBases());
        System.out.println("bestAlignment = " + bestAlignment);
        return 1;
    }


    /**
     * Provide an initial value for reduce computations.
     * @return Initial value of reduce.
     */
    @Override
    public Long reduceInit() {
        return 0L;
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     * @param value result of the map.
     * @param accum accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    @Override
    public Long reduce(Integer value, Long accum) {
        return value + accum;
    }

    @Override
    public void onTraversalDone(Long result) {
        aligner.close();
        super.onTraversalDone(result);
    }

}
