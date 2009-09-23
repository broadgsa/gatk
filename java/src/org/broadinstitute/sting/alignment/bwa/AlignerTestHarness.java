package org.broadinstitute.sting.alignment.bwa;

import org.broadinstitute.sting.alignment.bwa.bwt.*;
import org.broadinstitute.sting.alignment.Aligner;
import org.broadinstitute.sting.alignment.Alignment;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;

/**
 * A test harness to ensure that the perfect aligner works.
 *
 * @author mhanna
 * @version 0.1
 */
public class AlignerTestHarness {
    private static BWT bwt;
    private static SuffixArray suffixArray;

    public static void main( String argv[] ) throws FileNotFoundException {
        if( argv.length != 5 ) {
            System.out.println("PerfectAlignerTestHarness <fasta> <bwt> <rbwt> <sa> <bam>");
            System.exit(1);
        }

        File referenceFile = new File(argv[0]);
        File bwtFile = new File(argv[1]);
        File rbwtFile = new File(argv[2]);
        File reverseSuffixArrayFile = new File(argv[3]);
        File bamFile = new File(argv[4]);

        align(referenceFile,bwtFile,rbwtFile,reverseSuffixArrayFile,bamFile);
    }

    private static void align(File referenceFile, File bwtFile, File rbwtFile, File reverseSuffixArrayFile, File bamFile) throws FileNotFoundException {
        BWT bwt = new BWTReader(bwtFile).read();

        Aligner aligner = new BWAAligner(bwtFile,rbwtFile,reverseSuffixArrayFile);
        int count = 0;

        SAMFileReader reader = new SAMFileReader(bamFile);
        reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        for(SAMRecord read: reader) {
            count++;
            //if( count > 500 ) break;
            //if( count != 39 /*&& count != 15*/ ) continue;
            //if( !read.getReadName().endsWith("1507:1636#0") )
            //    continue;

            boolean skipRead = false;

            for( CigarElement cigarElement: read.getCigar().getCigarElements() ) {
                if( cigarElement.getOperator() != CigarOperator.M ) {
                    System.out.printf("Skipping read %s because it features indels%n", read.getReadName());
                    skipRead = true;
                }
            }

            if(read.getReadString().indexOf("N") >= 0) {
                System.out.printf("Skipping read %s because it contains Ns%n", read.getReadName());
                skipRead = true;
            }

            if(skipRead) continue;

            List<Alignment> alignments = aligner.align(read);
            if(alignments.size() == 0 )
                throw new StingException(String.format("Unable to align read %s to reference.",read.getReadName()));

            //System.out.printf("%s: Aligned read to reference at position %d with %d mismatches.%n", read.getReadName(), alignments.get(0).getAlignmentStart(), alignments.get(0).getScore());

            Alignment alignment = alignments.get(0);
            if( read.getAlignmentStart() != alignment.getAlignmentStart() ) {
                IndexedFastaSequenceFile reference = new IndexedFastaSequenceFile(referenceFile);
                String expectedRef = new String(reference.getSubsequenceAt(reference.getSequenceDictionary().getSequences().get(0).getSequenceName(),read.getAlignmentStart(),read.getAlignmentStart()+read.getReadLength()-1).getBases());
                int expectedMismatches = 0;
                for( int i = 0; i < read.getReadLength(); i++ ) {
                    if( read.getReadBases()[i] != expectedRef.charAt(i) )
                        expectedMismatches++;
                }

                String alignedRef = new String(reference.getSubsequenceAt(reference.getSequenceDictionary().getSequences().get(0).getSequenceName(),alignments.get(0).getAlignmentStart(),alignments.get(0).getAlignmentStart()+read.getReadLength()-1).getBases());
                int actualMismatches = 0;
                for( int i = 0; i < read.getReadLength(); i++ ) {
                    if( read.getReadBases()[i] != expectedRef.charAt(i) )
                        actualMismatches++;
                }

                if( expectedMismatches != actualMismatches ) {
                    System.out.printf("read          = %s%n", read.getReadString());
                    System.out.printf("expected ref  = %s%n", expectedRef);
                    System.out.printf("actual ref    = %s%n", alignedRef);
                    throw new StingException(String.format("Read %s was placed at incorrect location; target alignment = %d; actual alignment = %d%n",read.getReadName(),read.getAlignmentStart(),alignment.getAlignmentStart()));
                }
            }

            if( count % 1000 == 0 )
                System.out.printf("%d reads examined.%n",count);                
        }

        System.out.printf("%d reads examined.%n",count);
    }

}
