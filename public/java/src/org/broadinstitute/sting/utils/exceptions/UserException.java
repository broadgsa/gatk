/*
 * Copyright (c) 2010, The Broad Institute
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

package org.broadinstitute.sting.utils.exceptions;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;

/**
 * Represents the common user errors detected by Sting / GATK
 *
 * Root class for all GATK user errors, as well as the container for errors themselves
 *
 * User: depristo
 * Date: Sep 3, 2010
 * Time: 2:24:09 PM
 */
@DocumentedGATKFeature(
        groupName = "User exceptions",
        summary = "Exceptions caused by incorrect user behavior, such as bad files, bad arguments, etc." )
public class UserException extends ReviewedStingException {
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

    public static class MalformedGenomeLoc extends UserException {
        public MalformedGenomeLoc(String message, GenomeLoc loc) {
            super(String.format("Badly formed genome loc: %s: %s", message, loc));
        }

        public MalformedGenomeLoc(String message) {
            super(String.format("Badly formed genome loc: %s", message));
        }
    }

    public static class BadInput extends UserException {
        public BadInput(String message) {
            super(String.format("Bad input: %s", message));
        }
    }

    public static class NotSupportedInGATKLite extends UserException {
        public NotSupportedInGATKLite(String message) {
            super(String.format("GATK Lite does support all of the features of the full version: %s", message));
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
            super(String.format("Unknown tribble type %s: %s", type, message));
        }
    }


    public static class BadTmpDir extends UserException {
        public BadTmpDir(String message) {
            super(String.format("Failure working with the tmp directory %s. Override with -Djava.io.tmpdir=X on the command line to a bigger/better file system.  Exact error was %s", System.getProperties().get("java.io.tmpdir"), message));
        }
    }

    public static class TooManyOpenFiles extends UserException {
        public TooManyOpenFiles() {
            super(String.format("There was a failure because there are too many files open concurrently; your system's open file handle limit is too small.  See the unix ulimit command to adjust this limit"));
        }
    }

    public static class NotEnoughMemory extends UserException {
        public NotEnoughMemory() {
            super(String.format("There was a failure because you did not provide enough memory to run this program.  See the -Xmx JVM argument to adjust the maximum heap size provided to Java"));
        }
    }

    public static class ErrorWritingBamFile extends UserException {
        public ErrorWritingBamFile(String message) {
            super(String.format("An error occurred when trying to write the BAM file.  Usually this happens when there is not enough space in the directory to which the data is being written (generally the temp directory) or when your system's open file handle limit is too small.  To tell Java to use a bigger/better file system use -Djava.io.tmpdir=X on the command line.  The exact error was %s", message));
        }
    }

    public static class CouldNotReadInputFile extends UserException {
        public CouldNotReadInputFile(String message, Exception e) {
            super(String.format("Couldn't read file because %s caused by %s", message, getMessage(e)));
        }

        public CouldNotReadInputFile(File file) {
            super(String.format("Couldn't read file %s", file.getAbsolutePath()));
        }

        public CouldNotReadInputFile(File file, String message) {
            super(String.format("Couldn't read file %s because %s", file.getAbsolutePath(), message));
        }

        public CouldNotReadInputFile(File file, String message, Exception e) {
            super(String.format("Couldn't read file %s because %s with exception %s", file.getAbsolutePath(), message, getMessage(e)));
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
            super(String.format("Couldn't write file %s because %s with exception %s", file.getAbsolutePath(), message, getMessage(e)));
        }

        public CouldNotCreateOutputFile(File file, String message) {
            super(String.format("Couldn't write file %s because %s", file.getAbsolutePath(), message));
        }

        public CouldNotCreateOutputFile(String filename, String message, Exception e) {
            super(String.format("Couldn't write file %s because %s with exception %s", filename, message, getMessage(e)));
        }

        public CouldNotCreateOutputFile(File file, Exception e) {
            super(String.format("Couldn't write file %s because exception %s", file.getAbsolutePath(), getMessage(e)));
        }

        public CouldNotCreateOutputFile(String message, Exception e) {
            super(message, e);
        }
    }

    public static class MissortedBAM extends UserException {
        public MissortedBAM(SAMFileHeader.SortOrder order, File file, SAMFileHeader header) {
            super(String.format("Missorted Input SAM/BAM files: %s is must be sorted in %s order but order was: %s", file, order, header.getSortOrder()));
        }

        public MissortedBAM(SAMFileHeader.SortOrder order, String message) {
            super(String.format("Missorted Input SAM/BAM files: files are not sorted in %s order; %s", order, message));
        }

        public MissortedBAM(SAMFileHeader.SortOrder order, SAMRecord read, String message) {
            super(String.format("Missorted Input SAM/BAM file %s: file sorted in %s order but %s is required; %s",
                    read.getFileSource().getReader(), read.getHeader().getSortOrder(), order, message));
        }

        public MissortedBAM(String message) {
            super(String.format("Missorted Input SAM/BAM files: %s", message));
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
            super(String.format("SAM/BAM file %s is malformed: %s", source, message));
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
        public ReadMissingReadGroup(SAMRecord read) {
            super(read, String.format("Read %s is either missing the read group or its read group is not defined in the BAM header, both of which are required by the GATK.  Please use http://www.broadinstitute.org/gsa/wiki/index.php/ReplaceReadGroups to fix this problem", read.getReadName()));
        }
    }

    public static class VariantContextMissingRequiredField extends UserException {
        public VariantContextMissingRequiredField(String field, VariantContext vc) {
            super(String.format("Variant at %s:%d is is missing the required field %s", vc.getChr(), vc.getStart(), field));
        }
    }

    public static class MissortedFile extends UserException {
        public MissortedFile(File file, String message, Exception e) {
            super(String.format("Missorted Input file: %s is must be sorted in coordinate order. %s and got error %s", file, message, getMessage(e)));
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
            super(String.format("Input files %s and %s have incompatible contigs: %s.\n  %s contigs = %s\n  %s contigs = %s",
                    name1, name2, message, name1, ReadUtils.prettyPrintSequenceRecords(dict1), name2, ReadUtils.prettyPrintSequenceRecords(dict2)));
        }
    }

    public static class LexicographicallySortedSequenceDictionary extends UserException {
        public LexicographicallySortedSequenceDictionary(String name, SAMSequenceDictionary dict) {
            super(String.format("Lexicographically sorted human genome sequence detected in %s."
                    + "\nFor safety's sake the GATK requires human contigs in karyotypic order: 1, 2, ..., 10, 11, ..., 20, 21, 22, X, Y with M either leading or trailing these contigs."
                    + "\nThis is because all distributed GATK resources are sorted in karyotypic order, and your processing will fail when you need to use these files."
                    + "\nYou can use the ReorderSam utility to fix this problem: http://www.broadinstitute.org/gsa/wiki/index.php/ReorderSam"
                    + "\n  %s contigs = %s",
                    name, name, ReadUtils.prettyPrintSequenceRecords(dict)));
        }
    }



    public static class MissingWalker extends UserException {
        public MissingWalker(String walkerName, String message) {
            super(String.format("Walker %s is not available: %s", walkerName, message));
        }
    }

    public static class CannotExecuteQScript extends UserException {
        public CannotExecuteQScript(String message, Exception e) {
            super(String.format("Unable to execute QScript: " + message), e);
        }
    }

    public static class CouldNotCreateReferenceIndexFile extends UserException {
        public CouldNotCreateReferenceIndexFile(File f, Exception e) {
            this(f, "", e);
        }

        public CouldNotCreateReferenceIndexFile(File f, String message, Exception e) {
            super(String.format("Index file %s does not exist but could not be created because: %s. ", f, message)
                    + (e == null ? "" : getMessage(e)));
        }
    }

    public static class CouldNotCreateReferenceIndexFileBecauseOfLock extends UserException.CouldNotCreateReferenceIndexFile {
        public CouldNotCreateReferenceIndexFileBecauseOfLock(File f) {
            super(f, "could not be written because an exclusive file lock could not be obtained. " +
                    "If you are running multiple instances of GATK, another GATK process is " +
                    "probably creating this file now, and has locked it. Please wait until this process finishes " +
                    "and try again.", null);
        }
    }

    public static class UnreadableKeyException extends UserException {
        public UnreadableKeyException ( File f, Exception e ) {
            super(String.format("Key file %s cannot be read (possibly the key file is corrupt?). Error was: %s. " +
                                "Please see http://www.broadinstitute.org/gsa/wiki/index.php/Phone_home for help.",
                                f.getAbsolutePath(), getMessage(e)));
        }

        public UnreadableKeyException ( String message, Exception e ) {
            this(String.format("%s. Error was: %s", message, getMessage(e)));
        }

        public UnreadableKeyException ( String message ) {
            super(String.format("Key file cannot be read (possibly the key file is corrupt?): %s. " +
                                "Please see http://www.broadinstitute.org/gsa/wiki/index.php/Phone_home for help.",
                                message));
        }
    }

    public static class KeySignatureVerificationException extends UserException {
        public KeySignatureVerificationException ( File f ) {
            super(String.format("The signature in key file %s failed cryptographic verification. " +
                                "If this key was valid in the past, it's likely been revoked. " +
                                "Please see http://www.broadinstitute.org/gsa/wiki/index.php/Phone_home " +
                                "for help.",
                                f.getAbsolutePath()));
        }
    }
}
