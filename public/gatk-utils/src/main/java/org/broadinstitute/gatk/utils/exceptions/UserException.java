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

package org.broadinstitute.gatk.utils.exceptions;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.sam.ReadUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFIndexType;
import htsjdk.variant.variantcontext.VariantContext;

import java.io.File;

/**
 * Represents the common user errors detected by GATK
 *
 * Root class for all GATK user errors, as well as the container for errors themselves
 */
@DocumentedGATKFeature(
        groupName = HelpConstants.DOCS_CAT_USRERR,
        summary = "Errors caused by incorrect user behavior, such as bad files, bad arguments, etc." )
public class UserException extends ReviewedGATKException {
    /**
     * The URL where people can get help messages.  Printed when an error occurs
     */
    public static final String PHONE_HOME_DOCS_URL = "http://gatkforums.broadinstitute.org/discussion/1250/what-is-phone-home-and-how-does-it-affect-me#latest";

    public UserException(String msg) { super(msg); }
    public UserException(String msg, Throwable e) { super(msg, e); }
    private UserException(Throwable e) { super("", e); } // cannot be called, private access

    protected static String getMessage(Throwable t) {
        String message = t.getMessage();
        return message != null ? message : t.getClass().getName();
    }

    public static class CommandLineException extends UserException {
        public CommandLineException(String message) {
            super(String.format("Invalid command line: %s", message));
        }
    }

    public static class MalformedReadFilterException extends CommandLineException {
        public MalformedReadFilterException(String message) {
            super(String.format("Malformed read filter: %s",message));
        }
    }

    public static class IncompatibleReadFiltersException extends CommandLineException {
        public IncompatibleReadFiltersException(final String filter1, final String filter2) {
            super(String.format("Two read filters are enabled that are incompatible and cannot be used simultaneously: %s and %s", filter1, filter2));
        }
    }

    public static class MalformedWalkerArgumentsException extends CommandLineException {
        public MalformedWalkerArgumentsException(String message) {
            super(String.format("Malformed walker argument: %s",message));
        }
    }

    public static class UnsupportedCigarOperatorException extends UserException {
        public UnsupportedCigarOperatorException(final CigarOperator co, final SAMRecord read, final String message) {
            super(String.format(
                "Unsupported CIGAR operator %s in read %s at %s:%d. %s",
                co,
                read.getReadName(),
                read.getReferenceName(),
                read.getAlignmentStart(),
                message));
        }
    }


    public static class MalformedGenomeLoc extends UserException {
        public MalformedGenomeLoc(String message, GenomeLoc loc) {
            super(String.format("Badly formed genome location: %s: %s", message, loc));
        }

        public MalformedGenomeLoc(String message) {
            super(String.format("Badly formed genome location: %s", message));
        }
    }

    public static class BadInput extends UserException {
        public BadInput(String message) {
            super(String.format("Bad input: %s", message));
        }
    }

    // todo -- fix up exception cause passing
    public static class MissingArgument extends CommandLineException {
        public MissingArgument(String arg, String message) {
            super(String.format("Argument %s was missing: %s", arg, message));
        }
    }

    public static class BadArgumentValue extends CommandLineException {
        public BadArgumentValue(String arg, String message) {
            super(String.format("Argument %s has a bad value: %s", arg, message));
        }
    }

    public static class UnknownTribbleType extends CommandLineException {
        public UnknownTribbleType(String type, String message) {
            super(String.format("Unknown variant input file type %s: %s", type, message));
        }
    }


    public static class BadTmpDir extends UserException {
        public BadTmpDir(String message) {
            super(String.format("An error occurred while working with the tmp directory %s. You can specify -Djava.io.tmpdir=X on the command line (before the -jar argument) where X is a directory path, to use a more appropriate temporary directory. The exact error was %s", System.getProperties().get("java.io.tmpdir"), message));
        }
    }

    public static class TooManyOpenFiles extends UserException {
        public TooManyOpenFiles() {
            super(String.format("An error occurred because there were too many files open concurrently; your system's open file handle limit is probably too small.  See the unix ulimit command to adjust this limit or ask your system administrator for help."));
        }
    }

    public static class LocalParallelizationProblem extends UserException {
        public LocalParallelizationProblem(final File file) {
            super(String.format("An error occurred because temporary file %s could not be found while running the GATK with more than one thread. Possible causes for this problem include: your system's open file handle limit is too small, your output or temp directories do not have sufficient space, or your system experienced a temporary instability. Your system administrator can help you resolve these problems.", file.getAbsolutePath()));
        }
    }

    public static class NotEnoughMemory extends UserException {
        public NotEnoughMemory() {
            super(String.format("An error occurred because you did not provide enough memory to run this program. You can use the -Xmx argument (before the -jar argument) to adjust the maximum heap size provided to Java. Note that this is a JVM argument, not a GATK argument."));
        }
    }

    public static class ErrorWritingBamFile extends UserException {
        public ErrorWritingBamFile(String message) {
            super(String.format("An error occurred when trying to write the BAM file.  Usually this happens when there is not enough space in the directory to which the data is being written (generally the temp directory) or when your system's open file handle limit is too small.  Your system administrator can help you resolve these issues. If you know what temporary directory to use, you can specify it by adding -Djava.io.tmpdir=X to the command line (before the -jar argument), where X is the directory path. The exact error was %s", message));
        }
    }

    public static class NoSpaceOnDevice extends UserException {
        public NoSpaceOnDevice() {
            super("Writing failed because there is no space left on the disk or hard drive. Please make some space or specify a different location for writing output files.");
        }
    }

    public static class CouldNotReadInputFile extends UserException {
        public CouldNotReadInputFile(String message, Exception e) {
            super(String.format("Could not read file because %s caused by %s", message, getMessage(e)));
        }

        public CouldNotReadInputFile(File file) {
            super(String.format("Could not read file %s", file.getAbsolutePath()));
        }

        public CouldNotReadInputFile(File file, String message) {
            super(String.format("Could not read file %s because %s", file.getAbsolutePath(), message));
        }

        public CouldNotReadInputFile(String file, String message) {
            super(String.format("Could not read file %s because %s", file, message));
        }

        public CouldNotReadInputFile(File file, String message, Exception e) {
            super(String.format("Could not read file %s because %s with exception %s", file.getAbsolutePath(), message, getMessage(e)));
        }

        public CouldNotReadInputFile(File file, Exception e) {
            this(file, getMessage(e));
        }

        public CouldNotReadInputFile(String message) {
            super(message);
        }
    }


    public static class CouldNotCreateOutputFile extends UserException {
        public CouldNotCreateOutputFile(File file, String message, Exception e) {
            super(String.format("Could not write file %s because %s with exception %s", file.getAbsolutePath(), message, getMessage(e)));
        }

        public CouldNotCreateOutputFile(File file, String message) {
            super(String.format("Could not write file %s because %s", file.getAbsolutePath(), message));
        }

        public CouldNotCreateOutputFile(String filename, String message, Exception e) {
            super(String.format("Could not write file %s because %s with exception %s", filename, message, getMessage(e)));
        }

        public CouldNotCreateOutputFile(File file, Exception e) {
            super(String.format("Could not write file %s because exception %s", file.getAbsolutePath(), getMessage(e)));
        }

        public CouldNotCreateOutputFile(String message, Exception e) {
            super(message, e);
        }
    }

    public static class MissortedBAM extends UserException {
        public MissortedBAM(SAMFileHeader.SortOrder order, File file, SAMFileHeader header) {
            super(String.format("Missorted input SAM/BAM/CRAM files: %s must be sorted in %s order but order was: %s. Please see " + HelpConstants.forumPost("discussion/1317/collected-faqs-about-input-files-for-sequence-read-data-bam-cram") + "for more information.", file, order, header.getSortOrder()));
        }

        public MissortedBAM(SAMFileHeader.SortOrder order, String message) {
            super(String.format("Missorted input SAM/BAM/CRAM files: files are not sorted in %s order. Please see " + HelpConstants.forumPost("discussion/1317/collected-faqs-about-input-files-for-sequence-read-data-bam-cram") + "for more information. Error details: %s", order, message));
        }

        public MissortedBAM(SAMFileHeader.SortOrder order, SAMRecord read, String message) {
            super(String.format("Missorted input SAM/BAM/CRAM file %s: file sorted in %s order but %s is required. Please see " + HelpConstants.forumPost("discussion/1317/collected-faqs-about-input-files-for-sequence-read-data-bam-cram") + "for more information. Error details: %s",
                    read.getFileSource().getReader(), read.getHeader().getSortOrder(), order, message));
        }

        public MissortedBAM(String message) {
            super(String.format("Missorted input SAM/BAM/CRAM files. Please see " + HelpConstants.forumPost("discussion/1317/collected-faqs-about-input-files-for-sequence-read-data-bam-cram") + "for more information. Error details: %s", message));
        }
    }

    public static class MalformedBAM extends UserException {
        public MalformedBAM(SAMRecord read, String message) {
            this(read.getFileSource() != null ? read.getFileSource().getReader().toString() : "(none)", message);
        }

        public MalformedBAM(File file, String message) {
            this(file.toString(), message);
        }

        public MalformedBAM(String source, String message) {
            super(String.format("SAM/BAM/CRAM file %s is malformed. Please see " + HelpConstants.forumPost("discussion/1317/collected-faqs-about-input-files-for-sequence-read-data-bam-cram") + "for more information. Error details: %s", source, message));
        }
    }

    public static class MisencodedBAM extends UserException {
        public MisencodedBAM(SAMRecord read, String message) {
            this(read.getFileSource() != null ? read.getFileSource().getReader().toString() : "(none)", message);
        }

        public MisencodedBAM(String source, String message) {
            super(String.format("SAM/BAM/CRAM file %s appears to be using the wrong encoding for quality scores: %s. Please see https://www.broadinstitute.org/gatk/guide?id=6470 for more details and options related to this error.", source, message));
        }
    }

    public static class MalformedVCF extends UserException {
        public MalformedVCF(String message, String line) {
            super(String.format("The provided VCF file is malformed at line %s: %s", line, message));
        }

        public MalformedVCF(String message) {
            super(String.format("The provided VCF file is malformed: %s", message));
        }

        public MalformedVCF(String message, int lineNo) {
            super(String.format("The provided VCF file is malformed at approximately line number %d: %s", lineNo, message));
        }
    }

    public static class MalformedBCF2 extends UserException {
        public MalformedBCF2( String message ) {
            super(String.format("Malformed BCF2 file: %s", message));
        }
    }

    public static class MalformedVCFHeader extends UserException {
        public MalformedVCFHeader(String message) {
            super(String.format("The provided VCF file has a malformed header: %s", message));
        }
    }

    public static class ReadMissingReadGroup extends MalformedBAM {
        public ReadMissingReadGroup(final SAMRecord read) {
            super(read, String.format("Read %s is missing the read group (RG) tag, which is required by the GATK. Please see " + HelpConstants.forumPost("discussion/59/companion-utilities-replacereadgroups to fix this problem"), read.getReadName()));
        }
    }

    public static class ReadHasUndefinedReadGroup extends MalformedBAM {
        public ReadHasUndefinedReadGroup(final SAMRecord read, final String rgID) {
            super(read, String.format("Read %s uses a read group (%s) that is not defined in the BAM header, which is not valid.  Please see " + HelpConstants.forumPost("discussion/59/companion-utilities-replacereadgroups to fix this problem"), read.getReadName(), rgID));
        }
    }

    public static class VariantContextMissingRequiredField extends UserException {
        public VariantContextMissingRequiredField(String field, VariantContext vc) {
            super(String.format("Variant at %s:%d is is missing the required field %s.", vc.getChr(), vc.getStart(), field));
        }
    }

    public static class MissortedFile extends UserException {
        public MissortedFile(File file, String message, Exception e) {
            super(String.format("Missorted input file: %s is must be sorted in coordinate order. Please see " + HelpConstants.forumPost("discussion/1317/collected-faqs-about-input-files-for-sequence-read-data-bam-cram") + "for more information. Error details: %s and got error %s", file, message, getMessage(e)));
        }
    }

    public static class FailsStrictValidation extends UserException {
        public FailsStrictValidation(File f, String message) {
            super(String.format("File %s fails strict validation: %s", f.getAbsolutePath(), message));
        }
    }

    public static class MalformedFile extends UserException {
        public MalformedFile(String message) {
            super(String.format("Unknown file is malformed: %s", message));
        }

        public MalformedFile(String message, Exception e) {
            super(String.format("Unknown file is malformed: %s caused by %s", message, getMessage(e)));
        }

        public MalformedFile(File f, String message) {
            super(String.format("File %s is malformed: %s", f.getAbsolutePath(), message));
        }

        public MalformedFile(File f, String message, Exception e) {
            super(String.format("File %s is malformed: %s caused by %s", f.getAbsolutePath(), message, getMessage(e)));
        }

        public MalformedFile(String name, String message) {
            super(String.format("File associated with name %s is malformed: %s", name, message));
        }

        public MalformedFile(String name, String message, Exception e) {
            super(String.format("File associated with name %s is malformed: %s caused by %s", name, message, getMessage(e)));
        }
     }

    public static class CannotExecuteRScript extends UserException {
        public CannotExecuteRScript(String message) {
            super(String.format("Unable to execute RScript command: " + message));
        }
        public CannotExecuteRScript(String message, Exception e) {
            super(String.format("Unable to execute RScript command: " + message), e);
        }
    }

    public static class DeprecatedArgument extends CommandLineException {
        public DeprecatedArgument(String param, String doc) {
            super(String.format("The parameter %s is deprecated.  %s",param,doc));
        }
    }


    public static class IncompatibleSequenceDictionaries extends UserException {
        public IncompatibleSequenceDictionaries(String message, String name1, SAMSequenceDictionary dict1, String name2, SAMSequenceDictionary dict2) {
            super(String.format("Input files %s and %s have incompatible contigs. Please see " + HelpConstants.forumPost("discussion/63/input-files-have-incompatible-contigs") + "for more information. Error details: %s.\n  %s contigs = %s\n  %s contigs = %s",
                    name1, name2, message, name1, ReadUtils.prettyPrintSequenceRecords(dict1), name2, ReadUtils.prettyPrintSequenceRecords(dict2)));
        }
    }

    public static class LexicographicallySortedSequenceDictionary extends UserException {
        public LexicographicallySortedSequenceDictionary(String name, SAMSequenceDictionary dict) {
            super(String.format("Lexicographically sorted human genome sequence detected in %s. Please see " + HelpConstants.forumPost("discussion/58/companion-utilities-reordersam") + "for more information. Error details: %s contigs = %s",
                    name, name, ReadUtils.prettyPrintSequenceRecords(dict)));
        }
    }

    public static class DeprecatedWalker extends UserException {
        public DeprecatedWalker(String walkerName, String version) {
            super(String.format("Walker %s is no longer available in the GATK; it has been deprecated since version %s", walkerName, version));
        }
    }

    public static class DeprecatedAnnotation extends UserException {
        public DeprecatedAnnotation(String annotationName, String version) {
            super(String.format("Annotation %s is no longer available in the GATK; it has been deprecated since version %s", annotationName, version));
        }
    }

    public static class CannotExecuteQScript extends UserException {
        public CannotExecuteQScript(String message) {
            super(String.format("Unable to execute QScript: " + message));
        }
        public CannotExecuteQScript(String message, Exception e) {
            super(String.format("Unable to execute QScript: " + message), e);
        }
    }

    public static class CannotHandleGzippedRef extends UserException {
        public CannotHandleGzippedRef() {
            super("The GATK cannot process compressed (.gz) reference sequences. Please unzip the file and try again.  Sorry for the inconvenience.");
        }
    }

    public static class MissingReferenceFaiFile extends UserException {
        public MissingReferenceFaiFile( final File indexFile, final File fastaFile ) {
            super(String.format("Fasta index file %s for reference %s does not exist. Please see %s for help creating it.",
                                indexFile.getAbsolutePath(), fastaFile.getAbsolutePath(),
                                HelpConstants.forumPost("discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference")));
        }
    }

    public static class MissingReferenceDictFile extends UserException {
        public MissingReferenceDictFile( final File dictFile, final File fastaFile ) {
            super(String.format("Fasta dict file %s for reference %s does not exist. Please see %s for help creating it.",
                                dictFile.getAbsolutePath(), fastaFile.getAbsolutePath(),
                                HelpConstants.forumPost("discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference")));
        }
    }

    public static class UnreadableKeyException extends UserException {
        public UnreadableKeyException ( File f, Exception e ) {
            super(String.format("Key file %s cannot be read (possibly the key file is corrupt?). Error was: %s. " +
                                "Please see %s for help.",
                                f.getAbsolutePath(), getMessage(e), PHONE_HOME_DOCS_URL));
        }

        public UnreadableKeyException ( String message, Exception e ) {
            this(String.format("%s. Error was: %s", message, getMessage(e)));
        }

        public UnreadableKeyException ( String message ) {
            super(String.format("Key file cannot be read (possibly the key file is corrupt?): %s. " +
                                "Please see %s for help.",
                                message, PHONE_HOME_DOCS_URL));
        }
    }

    public static class KeySignatureVerificationException extends UserException {
        public KeySignatureVerificationException ( File f ) {
            super(String.format("The signature in key file %s failed cryptographic verification. " +
                                "If this key was valid in the past, it's likely been revoked. " +
                                "Please see %s for help.",
                                f.getAbsolutePath(), PHONE_HOME_DOCS_URL));
        }
    }

    public static class GVCFIndexException extends UserException {
        public GVCFIndexException (GATKVCFIndexType indexType, int indexParameter) {
            super(String.format("GVCF output requires a specific indexing strategy.  Please re-run including the arguments " +
                    "-variant_index_type %s -variant_index_parameter %d.",
                    indexType, indexParameter));
        }
    }

    /**
     * A special exception that happens only in the case where
     * the filesystem, by design or configuration, is completely unable
     * to handle locking.  This exception will specifically NOT be thrown
     * in the case where the filesystem handles locking but is unable to
     * acquire a lock due to concurrency.
     */
    public static class FileSystemInabilityToLockException extends UserException {
        public FileSystemInabilityToLockException( String message ) {
            super(message);
        }

        public FileSystemInabilityToLockException( String message, Exception innerException ) {
            super(message,innerException);
        }
    }

    public static class IncompatibleRecalibrationTableParameters extends UserException {
        public IncompatibleRecalibrationTableParameters(String s) {
            super(s);
        }
    }

    /**
     * A trivial specialization of  UserException to mark that a hardware feature is not supported
     */
    public static class HardwareFeatureException extends UserException {
        public HardwareFeatureException(String message) {
            super(message);
        }
    }
}
