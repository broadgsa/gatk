/*
 * Copyright (c) 2011, The Broad Institute
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

package org.broadinstitute.sting.utils.text;

import org.broadinstitute.sting.commandline.ParsingEngine;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.commandline.Tags;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.refdata.tracks.FeatureManager;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * A collection of convenience methods for working with list files.
 */
public class ListFileUtils {
    /**
     * Lines starting with this String in .list files are considered comments.
     */
    public static final String LIST_FILE_COMMENT_START = "#";        

    /**
     * Unpack the bam files to be processed, given a list of files.  That list of files can
     * itself contain entries which are lists of other files to be read (note: you cannot have lists
     * of lists of lists). Lines in .list files containing only whitespace or which begin with
     * LIST_FILE_COMMENT_START are ignored.
     *
     * @param samFiles The sam files, in string format.
     * @return a flattened list of the bam files provided
     */
    public static List<SAMReaderID> unpackBAMFileList(final List<String> samFiles, final ParsingEngine parser) {
        List<SAMReaderID> unpackedReads = new ArrayList<SAMReaderID>();
        for( String inputFileName: samFiles ) {
            Tags inputFileNameTags = parser.getTags(inputFileName);
            inputFileName = expandFileName(inputFileName);
            if (inputFileName.toLowerCase().endsWith(".list") ) {
                try {
                    for ( String fileName : new XReadLines(new File(inputFileName), true) ) {
                        if ( fileName.length() > 0 && ! fileName.startsWith(LIST_FILE_COMMENT_START) ) {
                            unpackedReads.add(new SAMReaderID(fileName,parser.getTags(inputFileName)));
                        }
                    }
                }
                catch( FileNotFoundException ex ) {
                    throw new UserException.CouldNotReadInputFile(new File(inputFileName), "Unable to find file while unpacking reads", ex);
                }
            }
            else if(inputFileName.toLowerCase().endsWith(".bam")) {
                unpackedReads.add(new SAMReaderID(inputFileName,inputFileNameTags));
            }
            else if(inputFileName.endsWith("stdin")) {
                unpackedReads.add(new SAMReaderID(inputFileName,inputFileNameTags));
            }
            else {
                throw new UserException.CommandLineException(String.format("The GATK reads argument (-I, --input_file) supports only BAM files with the .bam extension and lists of BAM files " +
                        "with the .list extension, but the file %s has neither extension.  Please ensure that your BAM file or list " +
                        "of BAM files is in the correct format, update the extension, and try again.",inputFileName));
            }
        }
        return unpackedReads;
    }

    /**
     * Convert command-line argument representation of ROD bindings to something more easily understandable by the engine.
     * @param RODBindings a text equivale
     * @return a list of expanded, bound RODs.
     */
    @Deprecated
    public static Collection<RMDTriplet> unpackRODBindingsOldStyle(final Collection<String> RODBindings, final ParsingEngine parser) {
        // todo -- this is a strange home for this code.  Move into ROD system
        Collection<RMDTriplet> rodBindings = new ArrayList<RMDTriplet>();

        for (String fileName: RODBindings) {
            final Tags tags = parser.getTags(fileName);
            fileName = expandFileName(fileName);

            List<String> positionalTags = tags.getPositionalTags();
            if(positionalTags.size() != 2)
                throw new UserException("Invalid syntax for -B (reference-ordered data) input flag.  " +
                        "Please use the following syntax when providing reference-ordered " +
                        "data: -B:<name>,<type> <filename>.");
            // Assume that if tags are present, those tags are name and type.
            // Name is always first, followed by type.
            String name = positionalTags.get(0);
            String type = positionalTags.get(1);

            RMDTriplet.RMDStorageType storageType = null;
            if(tags.getValue("storage") != null)
                storageType = Enum.valueOf(RMDTriplet.RMDStorageType.class,tags.getValue("storage"));
            else if(fileName.toLowerCase().endsWith("stdin"))
                storageType = RMDTriplet.RMDStorageType.STREAM;
            else
                storageType = RMDTriplet.RMDStorageType.FILE;

            rodBindings.add(new RMDTriplet(name,type,fileName,storageType,tags));
        }

        return rodBindings;
    }

    /**
     * Convert command-line argument representation of ROD bindings to something more easily understandable by the engine.
     * @param RODBindings a text equivale
     * @return a list of expanded, bound RODs.
     */
    public static Collection<RMDTriplet> unpackRODBindings(final Collection<RodBinding> RODBindings, final ParsingEngine parser) {
        // todo -- this is a strange home for this code.  Move into ROD system
        Collection<RMDTriplet> rodBindings = new ArrayList<RMDTriplet>();
        FeatureManager builderForValidation = new FeatureManager();

        for (RodBinding rodBinding: RODBindings) {
            String argValue = rodBinding.getSource();
            String fileName = expandFileName(argValue);
            String name = rodBinding.getName();
            String type = rodBinding.getTribbleType();

            RMDTriplet.RMDStorageType storageType = null;
            if(rodBinding.getTags().getValue("storage") != null)
                storageType = Enum.valueOf(RMDTriplet.RMDStorageType.class,rodBinding.getTags().getValue("storage"));
            else if(fileName.toLowerCase().endsWith("stdin"))
                storageType = RMDTriplet.RMDStorageType.STREAM;
            else
                storageType = RMDTriplet.RMDStorageType.FILE;

            RMDTriplet triplet = new RMDTriplet(name,type,fileName,storageType,rodBinding.getTags());

            // validate triplet type
            FeatureManager.FeatureDescriptor descriptor = builderForValidation.getByTriplet(triplet);
            if ( descriptor == null )
                throw new UserException.UnknownTribbleType(rodBinding.getTribbleType(),
                        String.format("Field %s had provided type %s but there's no such Tribble type.  The compatible types are: %n%s",
                                rodBinding.getName(), rodBinding.getTribbleType(), builderForValidation.userFriendlyListOfAvailableFeatures(rodBinding.getType())));
            if ( ! rodBinding.getType().isAssignableFrom(descriptor.getFeatureClass()) )
                throw new UserException.BadArgumentValue(rodBinding.getName(),
                        String.format("Field %s expects Features of type %s, but the input file produces Features of type %s. The compatible types are: %n%s",
                                rodBinding.getName(), rodBinding.getType().getSimpleName(), descriptor.getSimpleFeatureName(),
                                builderForValidation.userFriendlyListOfAvailableFeatures(rodBinding.getType())));


            rodBindings.add(triplet);
        }

        return rodBindings;
    }

    /**
     * Expand any special characters that appear in the filename.  Right now, '-' is expanded to
     * '/dev/stdin' only, but in the future, special characters like '~' and '*' that are passed
     * directly to the command line in some circumstances could be expanded as well.  Be careful
     * when adding UNIX-isms.
     * @param argument the text appearing on the command-line.
     * @return An expanded string suitable for opening by Java/UNIX file handling utilities.
     */
    private static String expandFileName(String argument) {
        if(argument.trim().equals("-"))
            return "/dev/stdin";
        return argument;
    }
}
