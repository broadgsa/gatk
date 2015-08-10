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

package org.broadinstitute.gatk.engine.alignment.bwa.java;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.*;
import org.broadinstitute.gatk.engine.alignment.Aligner;
import org.broadinstitute.gatk.engine.alignment.Alignment;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.io.File;
import java.io.FileNotFoundException;

/**
 * A test harness to ensure that the perfect aligner works.
 *
 * @author mhanna
 * @version 0.1
 */
public class AlignerTestHarness {
    public static void main( String argv[] ) throws FileNotFoundException {
        if( argv.length != 6 ) {
            System.out.println("PerfectAlignerTestHarness <fasta> <bwt> <rbwt> <sa> <rsa> <bam>");
            System.exit(1);
        }

        File referenceFile = new File(argv[0]);
        File bwtFile = new File(argv[1]);
        File rbwtFile = new File(argv[2]);
        File suffixArrayFile = new File(argv[3]);
        File reverseSuffixArrayFile = new File(argv[4]);
        File bamFile = new File(argv[5]);

        align(referenceFile,bwtFile,rbwtFile,suffixArrayFile,reverseSuffixArrayFile,bamFile);
    }

    private static void align(File referenceFile, File bwtFile, File rbwtFile, File suffixArrayFile, File reverseSuffixArrayFile, File bamFile) throws FileNotFoundException {
        Aligner aligner = new BWAJavaAligner(bwtFile,rbwtFile,suffixArrayFile,reverseSuffixArrayFile);
        int count = 0;

        SAMFileReader reader = new SAMFileReader(bamFile);
        reader.setValidationStringency(ValidationStringency.SILENT);

        int mismatches = 0;
        int failures = 0;

        for(SAMRecord read: reader) {
            count++;
            if( count > 200000 ) break;
            //if( count < 366000 ) continue;
            //if( count > 2 ) break;
            //if( !read.getReadName().endsWith("SL-XBC:1:82:506:404#0") )
            //    continue;
            //if( !read.getReadName().endsWith("SL-XBC:1:36:30:1926#0") )
            //    continue;
            //if( !read.getReadName().endsWith("SL-XBC:1:60:1342:1340#0") )
            //    continue;

            SAMRecord alignmentCleaned = null;
            try {
                alignmentCleaned = (SAMRecord)read.clone();
            }
            catch( CloneNotSupportedException ex ) {
                throw new ReviewedGATKException("SAMRecord clone not supported", ex);
            }

            if( alignmentCleaned.getReadNegativeStrandFlag() )
                alignmentCleaned.setReadBases(BaseUtils.simpleReverseComplement(alignmentCleaned.getReadBases()));

            alignmentCleaned.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
            alignmentCleaned.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
            alignmentCleaned.setMappingQuality(SAMRecord.NO_MAPPING_QUALITY);
            alignmentCleaned.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);

            // Clear everything except flags pertaining to pairing and set 'unmapped' status to true.
            alignmentCleaned.setFlags(alignmentCleaned.getFlags() & 0x00A1 | 0x000C);

            Iterable<Alignment[]> alignments = aligner.getAllAlignments(alignmentCleaned.getReadBases());
            if(!alignments.iterator().hasNext() ) {
                //throw new GATKException(String.format("Unable to align read %s to reference; count = %d",read.getReadName(),count));
                System.out.printf("Unable to align read %s to reference; count = %d%n",read.getReadName(),count);
                failures++;
            }

            Alignment foundAlignment = null;
            for(Alignment[] alignmentsOfQuality: alignments) {
                for(Alignment alignment: alignmentsOfQuality) {
                    if( read.getReadNegativeStrandFlag() != alignment.isNegativeStrand() )
                        continue;
                    if( read.getAlignmentStart() != alignment.getAlignmentStart() )
                        continue;

                    foundAlignment = alignment;                    
                }
            }

            if( foundAlignment != null ) {
                //System.out.printf("%s: Aligned read to reference at position %d with %d mismatches, %d gap opens, and %d gap extensions.%n", read.getReadName(), foundAlignment.getAlignmentStart(), foundAlignment.getMismatches(), foundAlignment.getGapOpens(), foundAlignment.getGapExtensions());
            }
            else {
                System.out.printf("Error aligning read %s%n", read.getReadName());

                mismatches++;

                IndexedFastaSequenceFile reference = new IndexedFastaSequenceFile(referenceFile);

                System.out.printf("read          = %s, position = %d, negative strand = %b%n", formatBasesBasedOnCigar(read.getReadString(),read.getCigar(),CigarOperator.DELETION),
                                                                                               read.getAlignmentStart(),
                                                                                               read.getReadNegativeStrandFlag());
                int numDeletions = numDeletionsInCigar(read.getCigar());
                String expectedRef = new String(reference.getSubsequenceAt(reference.getSequenceDictionary().getSequences().get(0).getSequenceName(),read.getAlignmentStart(),read.getAlignmentStart()+read.getReadLength()+numDeletions-1).getBases());
                System.out.printf("expected ref  = %s%n", formatBasesBasedOnCigar(expectedRef,read.getCigar(),CigarOperator.INSERTION));

                for(Alignment[] alignmentsOfQuality: alignments) {
                    for(Alignment alignment: alignmentsOfQuality) {
                        System.out.println();

                        Cigar cigar = ((BWAAlignment)alignment).getCigar();

                        System.out.printf("read          = %s%n", formatBasesBasedOnCigar(read.getReadString(),cigar,CigarOperator.DELETION));

                        int deletionCount = ((BWAAlignment)alignment).getNumberOfBasesMatchingState(AlignmentState.DELETION);
                        String alignedRef = new String(reference.getSubsequenceAt(reference.getSequenceDictionary().getSequences().get(0).getSequenceName(),alignment.getAlignmentStart(),alignment.getAlignmentStart()+read.getReadLength()+deletionCount-1).getBases());
                        System.out.printf("actual ref    = %s, position = %d, negative strand = %b%n", formatBasesBasedOnCigar(alignedRef,cigar,CigarOperator.INSERTION),
                                alignment.getAlignmentStart(),
                                alignment.isNegativeStrand());
                    }
                }

                //throw new GATKException(String.format("Read %s was placed at incorrect location; count = %d%n",read.getReadName(),count));                
            }


            if( count % 1000 == 0 )
                System.out.printf("%d reads examined.%n",count);                
        }

        System.out.printf("%d reads examined; %d mismatches; %d failures.%n",count,mismatches,failures);
    }

    private static String formatBasesBasedOnCigar( String bases, Cigar cigar, CigarOperator toBlank ) {
        StringBuilder formatted = new StringBuilder();
        int readIndex = 0;
        for(CigarElement cigarElement: cigar.getCigarElements()) {
            if(cigarElement.getOperator() == toBlank) {
                int number = cigarElement.getLength();
                while( number-- > 0 ) formatted.append(' ');
            }
            else {
                int number = cigarElement.getLength();
                while( number-- > 0 ) formatted.append(bases.charAt(readIndex++));
            }
        }
        return formatted.toString();
    }

    private static int numDeletionsInCigar( Cigar cigar ) {
        int numDeletions = 0;
        for(CigarElement cigarElement: cigar.getCigarElements()) {
            if(cigarElement.getOperator() == CigarOperator.DELETION)
                numDeletions += cigarElement.getLength();
        }
        return numDeletions;
    }
}
