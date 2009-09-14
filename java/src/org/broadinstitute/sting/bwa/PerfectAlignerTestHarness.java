package org.broadinstitute.sting.bwa;

import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.StingException;

import java.io.File;
import java.io.FileNotFoundException;

import net.sf.picard.reference.ReferenceSequence;

/**
 * A test harness to ensure that the perfect aligner works.
 *
 * @author mhanna
 * @version 0.1
 */
public class PerfectAlignerTestHarness {
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
        if( argv.length != 3 ) {
            System.out.println("PerfectAlignerTestHarness <fasta> <bwt> <sa>");
            System.exit(1);
        }

        File referenceFile = new File(argv[0]);
        IndexedFastaSequenceFile reference = new IndexedFastaSequenceFile(referenceFile);

        File bwtFile = new File(argv[1]);
        BWTReader reader = new BWTReader(bwtFile);
        bwt = reader.read();

        File suffixArrayFile = new File(argv[2]);
        SuffixArrayReader suffixArrayReader = new SuffixArrayReader(suffixArrayFile);
        suffixArray = suffixArrayReader.read();

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

    public static int align( String read ) {
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
