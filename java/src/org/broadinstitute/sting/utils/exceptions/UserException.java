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
import net.sf.samtools.SAMSequenceRecord;

import java.io.File;
import java.util.Arrays;

/**
 * Represents the common user errors detected by Sting / GATK
 *
 * Root class for all GATK user errors, as well as the container for errors themselves
 *
 * User: depristo
 * Date: Sep 3, 2010
 * Time: 2:24:09 PM
 */
public class UserException extends ReviewedStingException {
    public UserException(String msg) { super(msg); }
    public UserException(String msg, Throwable e) { super(msg, e); }
    private UserException(Throwable e) { super("", e); } // cannot be called, private access

    public static class CommandLineException extends UserException {
        public CommandLineException(String message) {
            super(String.format("Invalid command line: %s", message));
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

    public static class BadTmpDir extends UserException {
        public BadTmpDir(String message) {
            super(String.format("Failure working with the tmp directory %s. Override with -Djava.io.tmpdir=X on the command line to a bigger/better file system.  Exact error was %s", System.getProperties().get("java.io.tmpdir"), message));
        }
    }

    public static class CouldNotReadInputFile extends UserException {
        public CouldNotReadInputFile(String message, Exception e) {
            super(String.format("Couldn't read file because %s caused by %s", message, e.getMessage()));
        }

        public CouldNotReadInputFile(File file, String message) {
            super(String.format("Couldn't read file %s because %s", file.getAbsolutePath(), message));
        }

        public CouldNotReadInputFile(File file, String message, Exception e) {
            super(String.format("Couldn't read file %s because %s with exception %s", file.getAbsolutePath(), message, e.getMessage()));
        }

        public CouldNotReadInputFile(File file, Exception e) {
            this(file, e.getMessage());
        }
    }


    public static class CouldNotCreateOutputFile extends UserException {
        public CouldNotCreateOutputFile(File file, String message, Exception e) {
            super(String.format("Couldn't write file %s because %s with exception %s", file.getAbsolutePath(), message, e.getMessage()));
        }

        public CouldNotCreateOutputFile(String filename, String message, Exception e) {
            super(String.format("Couldn't write file %s because %s with exception %s", filename, message, e.getMessage()));
        }

        public CouldNotCreateOutputFile(File file, Exception e) {
            super(String.format("Couldn't write file %s because exception %s", file.getAbsolutePath(), e.getMessage()));
        }
    }

    public static class MalformedBam extends UserException {
        public MalformedBam(SAMRecord read, String message) {
            super(String.format("SAM/BAM file %s is malformed: %s", read.getFileSource().getReader(), message));
        }
    }

    public static class ReadMissingReadGroup extends MalformedBam {
        public ReadMissingReadGroup(SAMRecord read) {
            super(read, String.format("Read %s is either missing the read group or its read group is not defined in the BAM header, both of which are required by the GATK", read.getReadName()));
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

    public static class MissortedFile extends UserException {
        public MissortedFile(File file, String message, Exception e) {
            super(String.format("Missorted Input file: %s is must be sorted in coordinate order. %s and got error %s", file, message, e.getMessage()));
        }
    }

    public static class MalformedFile extends UserException {
        public MalformedFile(String message) {
            super(String.format("Unknown file is malformed: %s", message));
        }

        public MalformedFile(String message, Exception e) {
            super(String.format("Unknown file is malformed: %s caused by %s", message, e.getMessage()));
        }

        public MalformedFile(File f, String message) {
            super(String.format("File %s is malformed: %s", f.getAbsolutePath(), message));
        }

        public MalformedFile(File f, String message, Exception e) {
            super(String.format("File %s is malformed: %s caused by %s", f.getAbsolutePath(), message, e.getMessage()));
        }

        public MalformedFile(String name, String message, Exception e) {
            super(String.format("File associated with name %s is malformed: %s caused by %s", name, message, e.getMessage()));
        }
     }

    public static class CannotExecuteRScript extends UserException {
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
                    name1, name2, message, name1, prettyPrintSequenceRecords(dict1), name2, prettyPrintSequenceRecords(dict2)));
        }
    }

    public static class LexicographicallySortedSequenceDictionary extends UserException {
        public LexicographicallySortedSequenceDictionary(String name, SAMSequenceDictionary dict) {
            super(String.format("Lexicographically sorted human genome sequence detected in %s."
                    + "\nFor safety's sake the GATK requires human contigs in karyotypic order: 1, 2, ..., 10, 11, ..., 20, 21, 22, X, Y with M either leading or trailing these contigs."
                    + "\nThis is because all distributed GATK resources are sorted in karyotypic order, and your processing will fail when you need to use these files"
                    + "\n  %s contigs = %s",
                    name, name, prettyPrintSequenceRecords(dict)));
        }
    }

    private static String prettyPrintSequenceRecords(SAMSequenceDictionary sequenceDictionary) {
        String[] sequenceRecordNames = new String[sequenceDictionary.size()];
        int sequenceRecordIndex = 0;
        for (SAMSequenceRecord sequenceRecord : sequenceDictionary.getSequences())
            sequenceRecordNames[sequenceRecordIndex++] = sequenceRecord.getSequenceName();
        return Arrays.deepToString(sequenceRecordNames);
    }

    public static class MissingWalker extends UserException {
        public MissingWalker(String walkerName, String message) {
            super(String.format("Walker %s is not available: %s", walkerName, message));
        }
    }
}
