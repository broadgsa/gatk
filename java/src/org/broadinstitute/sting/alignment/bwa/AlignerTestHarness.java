package org.broadinstitute.sting.alignment.bwa;

import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.alignment.bwa.bwt.SuffixArrayReader;
import org.broadinstitute.sting.alignment.bwa.bwt.*;
import org.broadinstitute.sting.alignment.Aligner;
import org.broadinstitute.sting.alignment.Alignment;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader;

/**
 * A test harness to ensure that the perfect aligner works.
 *
 * @author mhanna
 * @version 0.1
 */
public class AlignerTestHarness {
    public static final String[] sampleReads = { "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGA",
                                                 "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTAG",
                                                 "GCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGCGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCAGAT",
                                                 "GCTTTTCATTCTGACTGCAACGGGCAATATGTATCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAA",
                                                 "GCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAA",
                                                 "CTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAAC",
                                                 "TTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTACTGAACT",
                                                 "TTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAAGAGTGTCTGATAGCAGCTTCTGAAC",
                                                 "TTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACT",
                                                 "CATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTT"};

    private static BWT bwt;
    private static SuffixArray suffixArray;

    public static void main( String argv[] ) throws FileNotFoundException {
        if( argv.length != 4 ) {
            System.out.println("PerfectAlignerTestHarness <fasta> <bwt> <rbwt>");
            System.exit(1);
        }

        File referenceFile = new File(argv[0]);
        File bwtFile = new File(argv[1]);
        File rbwtFile = new File(argv[2]);
        File bamFile = new File(argv[3]);
        File suffixArrayFile = null;

        align(referenceFile,bwtFile,rbwtFile,bamFile);
    }

    private static void align(File referenceFile, File bwtFile, File rbwtFile, File bamFile) throws FileNotFoundException {
        BWT bwt = new BWTReader(bwtFile).read();

        Aligner aligner = new BWAAligner(bwtFile,rbwtFile);
        int count = 0;

        SAMFileReader reader = new SAMFileReader(bamFile);
        reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        for(SAMRecord read: reader) {
            count++;
            if( count > 15 ) break;
            //if( count != 1 && count != 15 ) continue;

            List<Alignment> alignments = aligner.align(read);
            if(alignments.size() > 0 )
                System.out.printf("%s: Aligned read to reference with %d mismatches.%n", read.getReadName(), alignments.get(0).getScore());
            else
                System.out.printf("%s: Failed to align read to reference.%n", read.getReadName());
        }
    }

    private static void alignPerfect(File referenceFile, File bwtFile, File suffixArrayFile) throws FileNotFoundException
    {
        IndexedFastaSequenceFile reference = new IndexedFastaSequenceFile(referenceFile);
        BWT bwt = new BWTReader(bwtFile).read();
        SuffixArray suffixArray = new SuffixArrayReader(suffixArrayFile).read();

        for( String read: sampleReads ) {
            int alignmentStart = align(read);
            if( alignmentStart < 0 ) {
                System.out.printf("Unable to align read %s%n",read);
                continue;
            }
            ReferenceSequence subsequence = reference.getSubsequenceAt(reference.getSequenceDictionary().getSequences().get(0).getSequenceName(),alignmentStart,alignmentStart+read.length()-1);
            for( int i = 0; i < subsequence.length(); i++) {
                if( subsequence.getBases()[i] != read.charAt(i) )
                    throw new StingException("Read is not an exact match!  Alignment has failed!");
            }
            System.out.printf("Read %s aligned to position %d%n", read, alignmentStart);
        }        
    }

    private static int align(String read) {
        int lowerBound = 0, upperBound = bwt.length();
        for( int i = read.length()-1; i >= 0; i-- ) {
            Base base = Base.fromASCII((byte)read.charAt(i));
            lowerBound = bwt.counts(base) + bwt.occurrences(base,lowerBound-1)+1;
            upperBound = bwt.counts(base) + bwt.occurrences(base,upperBound);
            if( lowerBound > upperBound ) return -1;
        }
        int alignmentStart = suffixArray.sequence[lowerBound]+1;
        //System.out.printf("Read = %s; final bounds: (%d->%d); suffix array = %d%n",read,lowerBound,upperBound,alignmentStart);
        return alignmentStart;
    }


}
