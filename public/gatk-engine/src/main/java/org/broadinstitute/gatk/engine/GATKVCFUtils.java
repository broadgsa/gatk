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

package org.broadinstitute.gatk.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.index.*;
import htsjdk.tribble.index.interval.IntervalTreeIndex;
import htsjdk.tribble.index.linear.LinearIndex;
import org.apache.log4j.Logger;
import htsjdk.tribble.Feature;
import htsjdk.tribble.index.interval.IntervalIndexCreator;
import htsjdk.tribble.index.linear.LinearIndexCreator;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.gatk.utils.commandline.ArgumentTypeDescriptor;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.engine.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.gatk.utils.collections.Pair;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.utils.variant.GATKVCFIndexType;

import java.io.*;
import java.lang.reflect.Field;
import java.util.*;


/**
 * A set of GATK-specific static utility methods for common operations on VCF files/records.
 */
public class GATKVCFUtils {

    /**
     * Constructor access disallowed...static utility methods only!
     */
    private GATKVCFUtils() { }

    public static final Logger logger = Logger.getLogger(GATKVCFUtils.class);
    public final static String GATK_COMMAND_LINE_KEY = "GATKCommandLine";

    public final static GATKVCFIndexType DEFAULT_INDEX_TYPE = GATKVCFIndexType.DYNAMIC_SEEK;  // by default, optimize for seek time.  All indices prior to Nov 2013 used this type.
    public final static Integer DEFAULT_INDEX_PARAMETER = -1;           // the default DYNAMIC_SEEK does not use a parameter
    // as determined experimentally Nov-Dec 2013
    public final static GATKVCFIndexType DEFAULT_GVCF_INDEX_TYPE = GATKVCFIndexType.LINEAR;
    public final static Integer DEFAULT_GVCF_INDEX_PARAMETER = 128000;

    // GVCF file extensions
    public final static String GVCF_EXT = "g.vcf";
    public final static String GVCF_GZ_EXT = "g.vcf.gz";

    // Message for using the deprecated --variant_index_type or --variant_index_parameter arguments.
    public final static String DEPRECATED_INDEX_ARGS_MSG = "Naming your output file using the .g.vcf extension will automatically set the appropriate values " +
            " for --variant_index_type and --variant_index_parameter";

    /**
     * Gets the appropriately formatted header for a VCF file describing this GATK run
     *
     * @param header the existing VCFHeader that we will be adding this command line argument header line to.  Existing
     *               command line argument header lines will be used to generate a unique header line key.
     * @param engine the GATK engine that holds the walker name, GATK version, and other information
     * @param argumentSources contains information on the argument values provided to the GATK for converting to a
     *                        command line string.  Should be provided from the data in the parsing engine.  Can be
     *                        empty in which case the command line will be the empty string.
     * @return VCF header line describing this run of the GATK.
     */
    public static VCFHeaderLine getCommandLineArgumentHeaderLine(final VCFHeader header, final GenomeAnalysisEngine engine, final Collection<Object> argumentSources) {
        if ( engine == null ) throw new IllegalArgumentException("engine cannot be null");
        if ( argumentSources == null ) throw new IllegalArgumentException("argumentSources cannot be null");

        final Map<String, String> attributes = new LinkedHashMap<>();
        attributes.put("ID", engine.getWalkerName());
        attributes.put("Version", CommandLineGATK.getVersionNumber());
        final Date date = new Date();
        attributes.put("Date", date.toString());
        attributes.put("Epoch", Long.toString(date.getTime()));
        attributes.put("CommandLineOptions", engine.createApproximateCommandLineArgumentString(argumentSources.toArray()));

        // in case the walker name contains space, remove any spaces
        String key = getCommandLineKey(header, engine.getWalkerName().replaceAll("\\s", ""));
        return new VCFSimpleHeaderLine(key, attributes);
    }

    // create a unique command line argument header line key.  This method will look for existing
    // keys using the same walker name and append a count after it to make it unique.
    private static String getCommandLineKey(final VCFHeader header, final String walkerName) {
        final Iterator<VCFHeaderLine> existingMetaDataIterator = header.getMetaDataInInputOrder().iterator();

        // the command line argument keys are in the format GATK_COMMAND_LINE_KEY.(walker name)
        final String searchKey = String.format("%s.%s", GATK_COMMAND_LINE_KEY, walkerName);

        int commandLineKeyCount = 0;
        VCFHeaderLine line;
        while ( existingMetaDataIterator.hasNext() ) {
            line = existingMetaDataIterator.next();
            // if we find another key that starts with the same text as the walker
            if ( line.getKey().startsWith(searchKey) )
                commandLineKeyCount++;
        }

        // if there are no existing keys with this same walker name, then just return the
        // GATK_COMMAND_LINE_KEY.(walker name) format
        if ( commandLineKeyCount == 0 )
            return searchKey;
        // otherwise append the count associated with this new command (existing + 1)
        else
            return String.format("%s.%d", searchKey, commandLineKeyCount+1);
    }

    public static <T extends Feature> Map<String, VCFHeader> getVCFHeadersFromRods(GenomeAnalysisEngine toolkit, List<RodBinding<T>> rodBindings) {
        // Collect the eval rod names
        final Set<String> names = new TreeSet<String>();
        for ( final RodBinding<T> evalRod : rodBindings )
            names.add(evalRod.getName());
        return getVCFHeadersFromRods(toolkit, names);
    }

    public static Map<String, VCFHeader> getVCFHeadersFromRods(GenomeAnalysisEngine toolkit) {
        return getVCFHeadersFromRods(toolkit, (Collection<String>)null);
    }

    public static Map<String, VCFHeader> getVCFHeadersFromRods(GenomeAnalysisEngine toolkit, Collection<String> rodNames) {
        Map<String, VCFHeader> data = new HashMap<String, VCFHeader>();

        // iterate to get all of the sample names
        List<ReferenceOrderedDataSource> dataSources = toolkit.getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            // ignore the rod if it's not in our list
            if ( rodNames != null && !rodNames.contains(source.getName()) )
                continue;

            if ( source.getHeader() != null && source.getHeader() instanceof VCFHeader )
                data.put(source.getName(), (VCFHeader)source.getHeader());
        }

        return data;
    }

    public static Map<String,VCFHeader> getVCFHeadersFromRodPrefix(GenomeAnalysisEngine toolkit,String prefix) {
        Map<String, VCFHeader> data = new HashMap<String, VCFHeader>();

        // iterate to get all of the sample names
        List<ReferenceOrderedDataSource> dataSources = toolkit.getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            // ignore the rod if lacks the prefix
            if ( ! source.getName().startsWith(prefix) )
                continue;

            if ( source.getHeader() != null && source.getHeader() instanceof VCFHeader )
                data.put(source.getName(), (VCFHeader)source.getHeader());
        }

        return data;
    }

    /**
     * Gets the header fields from all VCF rods input by the user
     *
     * @param toolkit    GATK engine
     *
     * @return a set of all fields
     */
    public static Set<VCFHeaderLine> getHeaderFields(GenomeAnalysisEngine toolkit) {
        return getHeaderFields(toolkit, null);
    }

    /**
     * Gets the header fields from all VCF rods input by the user
     *
     * @param toolkit    GATK engine
     * @param rodNames   names of rods to use, or null if we should use all possible ones
     *
     * @return a set of all fields
     */
    public static Set<VCFHeaderLine> getHeaderFields(GenomeAnalysisEngine toolkit, Collection<String> rodNames) {

        // keep a map of sample name to occurrences encountered
        TreeSet<VCFHeaderLine> fields = new TreeSet<VCFHeaderLine>();

        // iterate to get all of the sample names
        List<ReferenceOrderedDataSource> dataSources = toolkit.getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            // ignore the rod if it's not in our list
            if ( rodNames != null && !rodNames.contains(source.getName()) )
                continue;

            if ( source.getRecordType().equals(VariantContext.class)) {
                VCFHeader header = (VCFHeader)source.getHeader();
                if ( header != null )
                    fields.addAll(header.getMetaDataInSortedOrder());
            }
        }

        return fields;
    }

    /**
     * Add / replace the contig header lines in the VCFHeader with the information in the GATK engine
     *
     * @param header the header to update
     * @param engine the GATK engine containing command line arguments and the master sequence dictionary
     */
    public static VCFHeader withUpdatedContigs(final VCFHeader header, final GenomeAnalysisEngine engine) {
        return VCFUtils.withUpdatedContigs(header, engine.getArguments().referenceFile, engine.getMasterSequenceDictionary());
    }

    /**
     * Create and return an IndexCreator
     * @param type
     * @param parameter
     * @param outFile
     * @return
     */
    public static IndexCreator getIndexCreator(GATKVCFIndexType type, int parameter, File outFile) {
        return getIndexCreator(type, parameter, outFile, null);
    }

    /**
     * Create and return an IndexCreator
     * @param type
     * @param parameter
     * @param outFile
     * @param sequenceDictionary
     * @return
     */
    public static IndexCreator getIndexCreator(GATKVCFIndexType type, int parameter, File outFile, SAMSequenceDictionary sequenceDictionary) {
        if (ArgumentTypeDescriptor.isCompressed(outFile.toString())) {
            if (type != GATKVCFUtils.DEFAULT_INDEX_TYPE || parameter != GATKVCFUtils.DEFAULT_INDEX_PARAMETER)
                logger.warn("Creating Tabix index for " + outFile + ", ignoring user-specified index type and parameter");

            if (sequenceDictionary == null)
                return new TabixIndexCreator(TabixFormat.VCF);
            else
                return new TabixIndexCreator(sequenceDictionary, TabixFormat.VCF);
        }

        IndexCreator idxCreator;
        switch (type) {
            case DYNAMIC_SEEK: idxCreator = new DynamicIndexCreator(outFile, IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME); break;
            case DYNAMIC_SIZE: idxCreator = new DynamicIndexCreator(outFile, IndexFactory.IndexBalanceApproach.FOR_SIZE); break;
            case LINEAR: idxCreator = new LinearIndexCreator(outFile, parameter); break;
            case INTERVAL: idxCreator = new IntervalIndexCreator(outFile, parameter); break;
            default: throw new IllegalArgumentException("Unknown IndexCreator type: " + type);
        }

        return idxCreator;
    }

    /**
     * Read all of the VCF records from source into memory, returning the header and the VariantContexts
     *
     * SHOULD ONLY BE USED FOR UNIT/INTEGRATION TESTING PURPOSES!
     *
     * @param source the file to read, must be in VCF4 format
     * @return
     * @throws java.io.IOException
     */
    public static Pair<VCFHeader, List<VariantContext>> readVCF(final File source) throws IOException {
        // read in the features
        final List<VariantContext> vcs = new ArrayList<VariantContext>();
        final VCFCodec codec = new VCFCodec();
        PositionalBufferedStream pbs = new PositionalBufferedStream(new FileInputStream(source));
        final LineIterator vcfSource = codec.makeSourceFromStream(pbs);
        try {
            final VCFHeader vcfHeader = (VCFHeader) codec.readActualHeader(vcfSource);

            while (vcfSource.hasNext()) {
                final VariantContext vc = codec.decode(vcfSource);
                if ( vc != null )
                    vcs.add(vc);
            }

            return new Pair<VCFHeader, List<VariantContext>>(vcfHeader, vcs);
        } finally {
            codec.close(vcfSource);
        }
    }

    /**
     * Check if the two indices are equivalent
     *
     * @param thisIndex index
     * @param otherIndex index
     * @return true if indices are equivalent, false otherwise.
     */
    public static boolean equivalentAbstractIndices(AbstractIndex thisIndex, AbstractIndex otherIndex){
        return thisIndex.getVersion() == otherIndex.getVersion() &&
                thisIndex.getIndexedFile().equals(otherIndex.getIndexedFile()) &&
                thisIndex.getIndexedFileSize() == otherIndex.getIndexedFileSize() &&
                thisIndex.getIndexedFileMD5().equals(otherIndex.getIndexedFileMD5()) &&
                thisIndex.getFlags() == otherIndex.getFlags();
    }

    /**
     * Check if the two indices are equivalent for a chromosome
     *
     * @param thisIndex index
     * @param otherIndex index
     * @param chr chromosome
     * @return true if indices are equivalent, false otherwise.
     * @throws NoSuchFieldException if index does not exist for a chromosome
     * @throws IllegalAccessException if index does not exist for a chromosome
     */
    public static boolean equivalentLinearIndices(LinearIndex thisIndex, LinearIndex otherIndex, String chr) throws NoSuchFieldException, IllegalAccessException {
        htsjdk.tribble.index.linear.LinearIndex.ChrIndex thisChr = (htsjdk.tribble.index.linear.LinearIndex.ChrIndex)getChrIndex(thisIndex, chr);
        htsjdk.tribble.index.linear.LinearIndex.ChrIndex otherChr = (htsjdk.tribble.index.linear.LinearIndex.ChrIndex)getChrIndex(otherIndex, chr);

        return  thisChr.getName().equals(otherChr.getName()) &&
                //thisChr.getTotalSize() == otherChr.getTotalSize() &&      TODO: why does this differ?
                thisChr.getNFeatures() == otherChr.getNFeatures() &&
                thisChr.getNBlocks() == otherChr.getNBlocks();
    }

    /**
     * Check if the two interval indices are equivalent for a chromosome
     *
     * @param thisIndex interval index
     * @param otherIndex interval index
     * @param chr chromosome
     * @return true if indices are equivalent, false otherwise.
     * @throws NoSuchFieldException if index does not exist for a chromosome
     * @throws IllegalAccessException if index does not exist for a chromosome
     */
    public static boolean equivalentIntervalIndices(IntervalTreeIndex thisIndex, IntervalTreeIndex otherIndex, String chr) throws NoSuchFieldException, IllegalAccessException {
        htsjdk.tribble.index.interval.IntervalTreeIndex.ChrIndex thisChr = (htsjdk.tribble.index.interval.IntervalTreeIndex.ChrIndex)getChrIndex(thisIndex, chr);
        htsjdk.tribble.index.interval.IntervalTreeIndex.ChrIndex otherChr = (htsjdk.tribble.index.interval.IntervalTreeIndex.ChrIndex)getChrIndex(otherIndex, chr);

        // TODO: compare trees?
        return thisChr.getName().equals(otherChr.getName());
    }

    /**
     * Get index for a chromosome
     *
     * @param index index
     * @param chr chromosome
     * @return index for the chromosome
     * @throws NoSuchFieldException if index does not exist for a chromosome
     * @throws IllegalAccessException if index does not exist for a chromosome
     */
    public static ChrIndex getChrIndex(AbstractIndex index, String chr) throws NoSuchFieldException, IllegalAccessException {
        Field f = AbstractIndex.class.getDeclaredField("chrIndices");
        f.setAccessible(true);
        LinkedHashMap<String, ChrIndex> chrIndices = (LinkedHashMap<String, ChrIndex>) f.get(index);
        return chrIndices.get(chr);
    }

    /**
     * Make an IndexCreator
     *
     * @param variantIndexType variant indexing strategy
     * @param variantIndexParameter variant indexing parameter
     * @param outputFile output variant file
     * @param sequenceDictionary collection of SAM sequence records
     * @return IndexCreator
     */
    public static IndexCreator makeIndexCreator(final GATKVCFIndexType variantIndexType, final int variantIndexParameter, final File outputFile, final SAMSequenceDictionary sequenceDictionary) {
        /*
        * If using the index arguments, log a warning.
        * If the genotype file has the GCVF extension (.g.vcf), use the default GCVF indexing.
        * Otherwise, use the default index type and parameter.
        */
        GATKVCFIndexType indexType = DEFAULT_INDEX_TYPE;
        int indexParameter = DEFAULT_INDEX_PARAMETER;
        if (usingNonDefaultIndexingArguments(variantIndexType, variantIndexParameter)) {
            indexType = variantIndexType;
            indexParameter = variantIndexParameter;
            logger.warn(DEPRECATED_INDEX_ARGS_MSG);
        } else if (outputFile.getName().endsWith("."  + GVCF_EXT) || outputFile.getName().endsWith("."  + GVCF_GZ_EXT)) {
            indexType = DEFAULT_GVCF_INDEX_TYPE;
            indexParameter = DEFAULT_GVCF_INDEX_PARAMETER;
        }

        return getIndexCreator(indexType, indexParameter, outputFile, sequenceDictionary);
    }

    /**
     * Check if not using the default indexing arguments' values
     *
     * @param variantIndexType variant indexing strategy
     * @param variantIndexParameter variant indexing parameter
     * @return true if the index type or parameter are not the default values, false otherwise
     */
    public static boolean usingNonDefaultIndexingArguments(final GATKVCFIndexType variantIndexType, final int variantIndexParameter) {
        return variantIndexType != GATKVCFUtils.DEFAULT_INDEX_TYPE || variantIndexParameter != GATKVCFUtils.DEFAULT_INDEX_PARAMETER;
    }

    /**
     * Check if using the GCVF indexing arguments' values
     *
     * @param variantIndexType variant indexing strategy
     * @param variantIndexParameter variant indexing parameter
     * @return true if the index type and parameter are the default GVCF values, false otherwise
     */
    public static boolean usingGVCFIndexingArguments(final GATKVCFIndexType variantIndexType, final int variantIndexParameter) {
        return variantIndexType == GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE && variantIndexParameter == GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER;
    }
}