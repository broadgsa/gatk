/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils.variant;

import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.FeatureCodecHeader;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.*;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;

/**
 * A set of GATK-specific static utility methods for common operations on VCF files/records.
 */
public class GATKVCFUtils {

    /**
     * Constructor access disallowed...static utility methods only!
     */
    private GATKVCFUtils() { }

    public final static String GATK_COMMAND_LINE_KEY = "GATKCommandLine";

    /**
     * Gets the appropriately formatted header for a VCF file describing this GATK run
     *
     * @param engine the GATK engine that holds the walker name, GATK version, and other information
     * @param argumentSources contains information on the argument values provided to the GATK for converting to a
     *                        command line string.  Should be provided from the data in the parsing engine.  Can be
     *                        empty in which case the command line will be the empty string.
     * @return VCF header line describing this run of the GATK.
     */
    public static VCFHeaderLine getCommandLineArgumentHeaderLine(final GenomeAnalysisEngine engine, final Collection<Object> argumentSources) {
        if ( engine == null ) throw new IllegalArgumentException("engine cannot be null");
        if ( argumentSources == null ) throw new IllegalArgumentException("argumentSources cannot be null");

        final Map<String, String> attributes = new LinkedHashMap<>();
        attributes.put("ID", engine.getWalkerName());
        attributes.put("Version", CommandLineGATK.getVersionNumber());
        final Date date = new Date();
        attributes.put("Date", date.toString());
        attributes.put("Epoch", Long.toString(date.getTime()));
        attributes.put("CommandLineOptions", engine.createApproximateCommandLineArgumentString(argumentSources.toArray()));
        return new VCFSimpleHeaderLine(GATK_COMMAND_LINE_KEY, attributes, Collections.<String>emptyList());
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
     * Utility class to read all of the VC records from a file
     *
     * @param source
     * @param codec
     * @return
     * @throws IOException
     */
    public final static Pair<VCFHeader, VCIterable> readAllVCs( final File source, final FeatureCodec<VariantContext> codec ) throws IOException {
        // read in the features
        PositionalBufferedStream pbs = new PositionalBufferedStream(new FileInputStream(source));
        FeatureCodecHeader header = codec.readHeader(pbs);
        pbs.close();

        pbs = new PositionalBufferedStream(new FileInputStream(source));
        pbs.skip(header.getHeaderEnd());

        final VCFHeader vcfHeader = (VCFHeader)header.getHeaderValue();
        return new Pair<VCFHeader, VCIterable>(vcfHeader, new VCIterable(pbs, codec, vcfHeader));
    }

    public static class VCIterable implements Iterable<VariantContext>, Iterator<VariantContext> {
        final PositionalBufferedStream pbs;
        final FeatureCodec<VariantContext> codec;
        final VCFHeader header;

        private VCIterable(final PositionalBufferedStream pbs, final FeatureCodec<VariantContext> codec, final VCFHeader header) {
            this.pbs = pbs;
            this.codec = codec;
            this.header = header;
        }

        @Override
        public Iterator<VariantContext> iterator() {
            return this;
        }

        @Override
        public boolean hasNext() {
            try {
                return ! pbs.isDone();
            } catch ( IOException e ) {
                throw new RuntimeException(e);
            }
        }

        @Override
        public VariantContext next() {
            try {
                final VariantContext vc = codec.decode(pbs);
                return vc == null ? null : vc.fullyDecode(header, false);
            } catch ( IOException e ) {
                throw new RuntimeException(e);
            }
        }

        @Override
        public void remove() {
        }
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
        FeatureCodecHeader header = codec.readHeader(pbs);
        pbs.close();

        pbs = new PositionalBufferedStream(new FileInputStream(source));
        pbs.skip(header.getHeaderEnd());

        final VCFHeader vcfHeader = (VCFHeader)header.getHeaderValue();

        while ( ! pbs.isDone() ) {
            final VariantContext vc = codec.decode(pbs);
            if ( vc != null )
                vcs.add(vc);
        }

        return new Pair<VCFHeader, List<VariantContext>>(vcfHeader, vcs);
    }
}