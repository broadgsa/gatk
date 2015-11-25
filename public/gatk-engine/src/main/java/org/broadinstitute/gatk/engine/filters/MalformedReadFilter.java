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

package org.broadinstitute.gatk.engine.filters;

import htsjdk.samtools.*;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.ReadProperties;
import org.broadinstitute.gatk.utils.ValidationExclusion;
import org.broadinstitute.gatk.engine.datasources.reads.SAMDataSource;
import org.broadinstitute.gatk.utils.exceptions.UserException;

/**
 * Filter out malformed reads
 *
 * <p>This filter is applied automatically by all GATK tools in order to protect them from crashing on reads that are
 * grossly malformed. There are a few issues (such as the absence of sequence bases) that will cause the run to fail with an
 * error, but these cases can be preempted by setting flags that cause the problem reads to also be filtered.</p>
 *
 * <h3>Usage example</h3>
 *
 * <h4>Set the malformed read filter to filter out reads that have no sequence bases</h4>
 * <pre>
 *     java -jar GenomeAnalysisTk.jar \
 *         -T ToolName \
 *         -R reference.fasta \
 *         -I input.bam \
 *         -o output.file \
 *         -filterNoBases
 * </pre>
 *
 * <p>Note that the MalformedRead filter itself does not need to be specified in the command line because it is set
 * automatically.</p>
 *
 * @author mhanna
 * @version 0.1
 */
public class MalformedReadFilter extends ReadFilter {


    private static final String FILTER_READS_WITH_N_CIGAR_ARGUMENT_FULL_NAME = "filter_reads_with_N_cigar" ;

    private SAMFileHeader header;

    @Argument(fullName = FILTER_READS_WITH_N_CIGAR_ARGUMENT_FULL_NAME, shortName = "filterRNC", doc = "Filter out reads with CIGAR containing the N operator, instead of failing with an error", required = false)
    boolean filterReadsWithNCigar = false;


    @Argument(fullName = "filter_mismatching_base_and_quals", shortName = "filterMBQ", doc = "Filter out reads with mismatching numbers of bases and base qualities, instead of failing with an error", required = false)
    boolean filterMismatchingBaseAndQuals = false;

    @Argument(fullName = "filter_bases_not_stored", shortName = "filterNoBases", doc = "Filter out reads with no stored bases (i.e. '*' where the sequence should be), instead of failing with an error", required = false)
    boolean filterBasesNotStored = false;

    /**
     * Indicates the applicable validation exclusions
     */
    private boolean allowNCigars;

    @Override
    public void initialize(final GenomeAnalysisEngine engine) {
        header = engine.getSAMFileHeader();
        ValidationExclusion validationExclusions = null;
        final SAMDataSource rds = engine.getReadsDataSource();
        if (rds != null) {
          final ReadProperties rps = rds.getReadsInfo();
          if (rps != null) {
            validationExclusions = rps.getValidationExclusionList();
          }
        }
        if (validationExclusions == null) {
            allowNCigars = false;
        } else {
            allowNCigars = validationExclusions.contains(ValidationExclusion.TYPE.ALLOW_N_CIGAR_READS);
        }
    }

    public boolean filterOut(final SAMRecord read) {
        // slowly changing the behavior to blow up first and filtering out if a parameter is explicitly provided
        return  !checkInvalidAlignmentStart(read) ||
                !checkInvalidAlignmentEnd(read) ||
                !checkAlignmentDisagreesWithHeader(this.header,read) ||
                !checkHasReadGroup(read) ||
                !checkMismatchingBasesAndQuals(read, filterMismatchingBaseAndQuals) ||
                !checkCigarDisagreesWithAlignment(read) ||
                !checkSeqStored(read, filterBasesNotStored) ||
                !checkCigarIsSupported(read,filterReadsWithNCigar,allowNCigars);
    }

    private static boolean checkHasReadGroup(final SAMRecord read) {
        if ( read.getReadGroup() == null ) {
            // there are 2 possibilities: either the RG tag is missing or it is not defined in the header
            final String rgID = (String)read.getAttribute(SAMTagUtil.getSingleton().RG);
            if ( rgID == null )
                throw new UserException.ReadMissingReadGroup(read);
            throw new UserException.ReadHasUndefinedReadGroup(read, rgID);
        }
        return true;
    }

    /**
     * Check for the case in which the alignment start is inconsistent with the read unmapped flag.
     * @param read The read to validate.
     * @return true if read start is valid, false otherwise.
     */
    private static boolean checkInvalidAlignmentStart(final SAMRecord read ) {
        // read is not flagged as 'unmapped', but alignment start is NO_ALIGNMENT_START
        if( !read.getReadUnmappedFlag() && read.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START )
            return false;
        // Read is not flagged as 'unmapped', but alignment start is -1
        if( !read.getReadUnmappedFlag() && read.getAlignmentStart() == -1 )
            return false;
        return true;
    }

    /**
     * Check for invalid end of alignments.
     * @param read The read to validate.
     * @return true if read end is valid, false otherwise.
     */
    private static boolean checkInvalidAlignmentEnd(final SAMRecord read ) {
        // Alignment aligns to negative number of bases in the reference.
        if( !read.getReadUnmappedFlag() && read.getAlignmentEnd() != -1 && (read.getAlignmentEnd()-read.getAlignmentStart()+1)<0 )
            return false;
        return true;
    }

    /**
     * Check to ensure that the alignment makes sense based on the contents of the header.
     * @param header The SAM file header.
     * @param read The read to verify.
     * @return true if alignment agrees with header, false othrewise.
     */
    private static boolean checkAlignmentDisagreesWithHeader(final SAMFileHeader header, final SAMRecord read ) {
        // Read is aligned to nonexistent contig
        if( read.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && read.getAlignmentStart() != SAMRecord.NO_ALIGNMENT_START )
            return false;
        final SAMSequenceRecord contigHeader = header.getSequence( read.getReferenceIndex() );
        // Read is aligned to a point after the end of the contig
        if( !read.getReadUnmappedFlag() && read.getAlignmentStart() > contigHeader.getSequenceLength() )
            return false;
        return true;
    }

    /**
     * Check for inconsistencies between the cigar string and the
     * @param read The read to validate.
     * @return true if cigar agrees with alignment, false otherwise.
     */
    private static boolean checkCigarDisagreesWithAlignment(final SAMRecord read) {
        // Read has a valid alignment start, but the CIGAR string is empty
        if( !read.getReadUnmappedFlag() &&
            read.getAlignmentStart() != -1 &&
            read.getAlignmentStart() != SAMRecord.NO_ALIGNMENT_START &&
            read.getAlignmentBlocks().size() < 0 )
            return false;
        return true;
    }

    /**
     * Check for unsupported CIGAR operators.
     * Currently the N operator is not supported.
     * @param read The read to validate.
     * @param filterReadsWithNCigar whether the offending read should just
     *                              be silently filtered or not.
     * @param allowNCigars whether reads that contain N operators in their CIGARs
     *                     can be processed or an exception should be thrown instead.
     * @throws UserException.UnsupportedCigarOperatorException
     *   if {@link #filterReadsWithNCigar} is <code>false</code> and
     *   the input read has some unsupported operation.
     * @return <code>true</code> if the read CIGAR operations are
     * fully supported, otherwise <code>false</code>, as long as
     * no exception has been thrown.
     */
    private static boolean checkCigarIsSupported(final SAMRecord read, final boolean filterReadsWithNCigar, final boolean allowNCigars) {
        if( containsNOperator(read)) {
            if (! filterReadsWithNCigar && !allowNCigars) {
                throw new UserException.UnsupportedCigarOperatorException(
                        CigarOperator.N,read,
                        "Perhaps you are"
                        + " trying to use RNA-Seq data?"
                        + " While we are currently actively working to"
                        + " support this data type unfortunately the"
                        + " GATK cannot be used with this data in its"
                        + " current form. You have the option of either"
                        + " filtering out all reads with operator "
                        + CigarOperator.N + " in their CIGAR string"
                        + " (please add --"
                        +  FILTER_READS_WITH_N_CIGAR_ARGUMENT_FULL_NAME
                        + " to your command line) or"
                        + " assume the risk of processing those reads as they"
                        + " are including the pertinent unsafe flag (please add -U"
                        + ' ' + ValidationExclusion.TYPE.ALLOW_N_CIGAR_READS
                        + " to your command line). Notice however that if you were"
                        + " to choose the latter, an unspecified subset of the"
                        + " analytical outputs of an unspecified subset of the tools"
                        + " will become unpredictable. Consequently the GATK team"
                        + " might well not be able to provide you with the usual support"
                        + " with any issue regarding any output");
            }
            return ! filterReadsWithNCigar;
        }
        return true;
    }

    private static boolean containsNOperator(final SAMRecord read) {
        final Cigar cigar = read.getCigar();
        if (cigar == null)   {
            return false;
        }
        for (final CigarElement ce : cigar.getCigarElements()) {
            if (ce.getOperator() == CigarOperator.N) {
                return true;
            }
        }
        return false;
    }

    /**
     * Check if the read has the same number of bases and base qualities
     * @param read the read to validate
     * @return true if they have the same number. False otherwise.
     */
    private static boolean checkMismatchingBasesAndQuals(final SAMRecord read, final boolean filterMismatchingBaseAndQuals) {
        final boolean result;
        if (read.getReadLength() == read.getBaseQualities().length)
            result = true;
        else if (filterMismatchingBaseAndQuals)
            result = false;
        else
            throw new UserException.MalformedBAM(read,
                    String.format("BAM file has a read with mismatching number of bases and base qualities. Offender: %s [%d bases] [%d quals].%s",
                            read.getReadName(), read.getReadLength(), read.getBaseQualities().length,
                            read.getBaseQualities().length == 0 ? " You can use --defaultBaseQualities to assign a default base quality for all reads, but this can be dangerous in you don't know what you are doing." : ""));

        return result;
    }

    /**
     * Check if the read has its base sequence stored
     * @param read the read to validate
     * @return true if the sequence is stored and false otherwise ("*" in the SEQ field).
     */
    protected static boolean checkSeqStored(final SAMRecord read, final boolean filterBasesNotStored) {

        if ( read.getReadBases() != SAMRecord.NULL_SEQUENCE )
            return true;

        if ( filterBasesNotStored )
            return false;

        throw new UserException.MalformedBAM(read, String.format("the BAM file has a read with no stored bases (i.e. it uses '*') which is not supported in the GATK; see the --filter_bases_not_stored argument. Offender: %s", read.getReadName()));
    }
}
