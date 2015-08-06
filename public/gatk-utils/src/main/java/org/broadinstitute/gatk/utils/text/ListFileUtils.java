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

package org.broadinstitute.gatk.utils.text;

import org.broadinstitute.gatk.utils.commandline.ParsingEngine;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.commandline.Tags;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;
import org.broadinstitute.gatk.utils.refdata.tracks.FeatureManager;
import org.broadinstitute.gatk.utils.refdata.utils.RMDTriplet;
import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;
import java.util.regex.Pattern;

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
     * @param parser Parser
     * @return a flattened list of the bam files provided
     */
    public static List<SAMReaderID> unpackBAMFileList(final List<String> samFiles, final ParsingEngine parser) {
        List<SAMReaderID> unpackedReads = new ArrayList<SAMReaderID>();
        for( String inputFileName: samFiles ) {
            Tags inputFileNameTags = parser.getTags(inputFileName);
            inputFileName = expandFileName(inputFileName);
            if (inputFileName.toLowerCase().endsWith(".list") ) {
                try {
                    for ( String fileName : new XReadLines(new File(inputFileName), true, LIST_FILE_COMMENT_START) ) {
                        unpackedReads.add(new SAMReaderID(fileName,parser.getTags(inputFileName)));
                    }
                }
                catch( FileNotFoundException ex ) {
                    throw new UserException.CouldNotReadInputFile(new File(inputFileName), "Unable to find file while unpacking reads", ex);
                }
            }
            else if(inputFileName.toLowerCase().endsWith(".bam") || inputFileName.toLowerCase().endsWith(".cram")) {
                unpackedReads.add(new SAMReaderID(inputFileName,inputFileNameTags));
            }
            else if(inputFileName.endsWith("stdin")) {
                unpackedReads.add(new SAMReaderID(inputFileName,inputFileNameTags));
            }
            else {
                throw new UserException.CommandLineException(String.format("The GATK reads argument (-I, --input_file) supports only BAM/CRAM files with the .bam/.cram extension and lists of BAM/CRAM files " +
                        "with the .list extension, but the file %s has neither extension.  Please ensure that your BAM/CRAM file or list " +
                        "of BAM/CRAM files is in the correct format, update the extension, and try again.",inputFileName));
            }
        }
        return unpackedReads;
    }

    /**
     * Convert command-line argument representation of ROD bindings to something more easily understandable by the engine.
     * @param RODBindings a text equivale
     * @param parser Parser
     * @return a list of expanded, bound RODs.
     */
    @Deprecated
    @SuppressWarnings("unused") // TODO: Who is still using this? External walkers?
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

            RMDTriplet.RMDStorageType storageType;
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
     * @param parser Parser
     * @return a list of expanded, bound RODs.
     */
    @SuppressWarnings("unchecked")
    public static Collection<RMDTriplet> unpackRODBindings(final Collection<RodBinding> RODBindings, @SuppressWarnings("unused") final ParsingEngine parser) {
        // todo -- this is a strange home for this code.  Move into ROD system
        Collection<RMDTriplet> rodBindings = new ArrayList<RMDTriplet>();
        FeatureManager builderForValidation = new FeatureManager();

        for (RodBinding rodBinding: RODBindings) {
            String argValue = rodBinding.getSource();
            String fileName = expandFileName(argValue);
            String name = rodBinding.getName();
            String type = rodBinding.getTribbleType();

            RMDTriplet.RMDStorageType storageType;
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

    /**
     * Returns a new set of values, containing a final set of values expanded from values
     * <p/>
     * Each element E of values can either be a literal string or a file ending in .list.
     * For each E ending in .list we try to read a file named E from disk, and if possible
     * all lines from that file are expanded into unique values.
     *
     * @param values Original values
     * @return entries from values or the files listed in values
     */
    public static Set<String> unpackSet(Collection<String> values) {
        if (values == null)
            throw new NullPointerException("values cannot be null");
        Set<String> unpackedValues = new LinkedHashSet<String>();
        // Let's first go through the list and see if we were given any files.
        // We'll add every entry in the file to our set, and treat the entries as
        // if they had been specified on the command line.
        for (String value : values) {
            File file = new File(value);
            if (value.toLowerCase().endsWith(".list") && file.exists()) {
                try {
                    unpackedValues.addAll(new XReadLines(file, true, LIST_FILE_COMMENT_START).readLines());
                } catch (IOException e) {
                    throw new UserException.CouldNotReadInputFile(file, e);
                }
            } else {
                unpackedValues.add(value);
            }
        }
        return unpackedValues;
    }

    /**
     * Returns a new set of values including only values listed by filters
     * <p/>
     * Each element E of values can either be a literal string or a file.  For each E,
     * we try to read a file named E from disk, and if possible all lines from that file are expanded
     * into unique names.
     * <p/>
     * Filters may also be a file of filters.
     *
     * @param values     Values or files with values
     * @param filters    Filters or files with filters
     * @param exactMatch If true match filters exactly, otherwise use as both exact and regular expressions
     * @return entries from values or the files listed in values, filtered by filters
     */
    public static Set<String> includeMatching(Collection<String> values, Collection<String> filters, boolean exactMatch) {
        return includeMatching(values, IDENTITY_STRING_CONVERTER, filters, exactMatch);
    }

    /**
     * Converts a type T to a String representation.
     *
     * @param <T> Type to convert to a String.
     */
    public static interface StringConverter<T> {
        String convert(T value);
    }

    /**
     * Returns a new set of values including only values matching filters
     * <p/>
     * Filters may also be a file of filters.
     * <p/>
     * The converter should convert T to a unique String for each value in the set.
     *
     * @param values     Values or files with values
     * @param converter  Converts values to strings
     * @param filters    Filters or files with filters
     * @param exactMatch If true match filters exactly, otherwise use as both exact and regular expressions
     * @return entries from values including only values matching filters
     */
    public static <T> Set<T> includeMatching(Collection<T> values, StringConverter<T> converter, Collection<String> filters, boolean exactMatch) {
        if (values == null)
            throw new NullPointerException("values cannot be null");
        if (converter == null)
            throw new NullPointerException("converter cannot be null");
        if (filters == null)
            throw new NullPointerException("filters cannot be null");

        Set<String> unpackedFilters = unpackSet(filters);
        Set<T> filteredValues = new LinkedHashSet<T>();
        Collection<Pattern> patterns = null;
        if (!exactMatch)
            patterns = compilePatterns(unpackedFilters);
        for (T value : values) {
            String converted = converter.convert(value);
            if (unpackedFilters.contains(converted)) {
                filteredValues.add(value);
            } else if (!exactMatch) {
                for (Pattern pattern : patterns)
                    if (pattern.matcher(converted).find())
                        filteredValues.add(value);
            }
        }
        return filteredValues;
    }
    
    /**
     * Returns a new set of values excluding any values matching filters.
     * <p/>
     * Filters may also be a file of filters.
     * <p/>
     * The converter should convert T to a unique String for each value in the set.
     *
     * @param values     Values or files with values
     * @param converter  Converts values to strings
     * @param filters    Filters or files with filters
     * @param exactMatch If true match filters exactly, otherwise use as both exact and regular expressions
     * @return entries from values exluding any values matching filters
     */
    public static <T> Set<T> excludeMatching(Collection<T> values, StringConverter<T> converter, Collection<String> filters, boolean exactMatch) {
        if (values == null)
            throw new NullPointerException("values cannot be null");
        if (converter == null)
            throw new NullPointerException("converter cannot be null");
        if (filters == null)
            throw new NullPointerException("filters cannot be null");

        Set<String> unpackedFilters = unpackSet(filters);
        Set<T> filteredValues = new LinkedHashSet<T>();
        filteredValues.addAll(values);
        Collection<Pattern> patterns = null;
        if (!exactMatch)
            patterns = compilePatterns(unpackedFilters);
        for (T value : values) {
            String converted = converter.convert(value);
            if (unpackedFilters.contains(converted)) {
                filteredValues.remove(value);
            } else if (!exactMatch) {
                for (Pattern pattern : patterns)
                    if (pattern.matcher(converted).find())
                        filteredValues.remove(value);
            }
        }
        return filteredValues;
    }

    private static Collection<Pattern> compilePatterns(Collection<String> filters) {
        Collection<Pattern> patterns = new ArrayList<Pattern>();
        for (String filter: filters) {
            patterns.add(Pattern.compile(filter));
        }
        return patterns;
    }

    protected static final StringConverter<String> IDENTITY_STRING_CONVERTER = new StringConverter<String>() {
        @Override
        public String convert(String value) {
            return value;
        }
    };
}
