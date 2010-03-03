package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.alignment.bwa.BWAAligner;
import org.broadinstitute.sting.alignment.bwa.BWAConfiguration;
import org.broadinstitute.sting.alignment.bwa.c.BWACAligner;
import org.broadinstitute.sting.alignment.Alignment;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.StringUtil;
import net.sf.picard.reference.ReferenceSequence;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Scanner;
import java.util.SortedSet;
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
            throw new StingException("Unable to load indel calls.");
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
                throw new StingException("Invalid indel type: " + type);                

            aligners.put(GenomeLocParser.createGenomeLoc(contig,start,stop),new BWACAligner(revisedReference,new BWAConfiguration()));
            if(++numAlignersCreated % 100 == 0)
                out.printf("Created %d aligners in %dms%n",++numAlignersCreated,System.currentTimeMillis()-startTime);
        }
    }

    @Override
    public Integer map(char[] ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
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
