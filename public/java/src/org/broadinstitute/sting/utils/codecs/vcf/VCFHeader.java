/*
 * Copyright (c) 2012, The Broad Institute
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

package org.broadinstitute.sting.utils.codecs.vcf;

import org.broad.tribble.util.ParsingUtils;

import java.util.*;


/**
 * @author aaron
 *         <p/>
 *         Class VCFHeader
 *         <p/>
 *         A class representing the VCF header
 */
public class VCFHeader {

    // the mandatory header fields
    public enum HEADER_FIELDS {
        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
    }

    // the associated meta data
    private final Set<VCFHeaderLine> mMetaData = new TreeSet<VCFHeaderLine>();
    private final Map<String, VCFInfoHeaderLine> mInfoMetaData = new HashMap<String, VCFInfoHeaderLine>();
    private final Map<String, VCFFormatHeaderLine> mFormatMetaData = new HashMap<String, VCFFormatHeaderLine>();
    private final Map<String, VCFHeaderLine> mOtherMetaData = new HashMap<String, VCFHeaderLine>();
    private final List<VCFContigHeaderLine> contigMetaData = new ArrayList<VCFContigHeaderLine>();

    // the list of auxillary tags
    private final Set<String> mGenotypeSampleNames = new LinkedHashSet<String>();

    // the character string that indicates meta data
    public static final String METADATA_INDICATOR = "##";

    // the header string indicator
    public static final String HEADER_INDICATOR = "#";

    public static final String SOURCE_KEY = "source";
    public static final String REFERENCE_KEY = "reference";
    public static final String CONTIG_KEY = "contig";
    public static final String INTERVALS_KEY = "intervals";

    // were the input samples sorted originally (or are we sorting them)?
    private boolean samplesWereAlreadySorted = true;

    // cache for efficient conversion of VCF -> VariantContext
    protected ArrayList<String> sampleNamesInOrder = null;
    protected HashMap<String, Integer> sampleNameToOffset = null;

    private boolean writeEngineHeaders = true;
    private boolean writeCommandLine = true;

    /**
     * create a VCF header, given a list of meta data and auxillary tags
     *
     * @param metaData     the meta data associated with this header
     */
    public VCFHeader(Set<VCFHeaderLine> metaData) {
        mMetaData.addAll(metaData);
        loadVCFVersion();
        loadMetaDataMaps();
    }

    /**
     * create a VCF header, given a list of meta data and auxillary tags
     *
     * @param metaData            the meta data associated with this header
     * @param genotypeSampleNames the sample names
     */
    public VCFHeader(Set<VCFHeaderLine> metaData, Set<String> genotypeSampleNames) {
        this(metaData);
        mGenotypeSampleNames.addAll(genotypeSampleNames);
        samplesWereAlreadySorted = ParsingUtils.isSorted(genotypeSampleNames);
        buildVCFReaderMaps(genotypeSampleNames);
    }

    /**
     * Tell this VCF header to use pre-calculated sample name ordering and the
     * sample name -> offset map.  This assumes that all VariantContext created
     * using this header (i.e., read by the VCFCodec) will have genotypes
     * occurring in the same order
     *
     * @param genotypeSampleNamesInAppearenceOrder genotype sample names, must iterator in order of appearence
     */
    private void buildVCFReaderMaps(Collection<String> genotypeSampleNamesInAppearenceOrder) {
        sampleNamesInOrder = new ArrayList<String>(genotypeSampleNamesInAppearenceOrder.size());
        sampleNameToOffset = new HashMap<String, Integer>(genotypeSampleNamesInAppearenceOrder.size());

        int i = 0;
        for ( final String name : genotypeSampleNamesInAppearenceOrder ) {
            sampleNamesInOrder.add(name);
            sampleNameToOffset.put(name, i++);
        }
        Collections.sort(sampleNamesInOrder);
    }


    /**
     * Adds a header line to the header metadata.
     *
     * @param headerLine Line to add to the existing metadata component.
     */
    public void addMetaDataLine(VCFHeaderLine headerLine) {
        mMetaData.add(headerLine);
    }

    /**
     * @return all of the VCF header lines of the ##contig form in order, or an empty set if none were present
     */
    public List<VCFContigHeaderLine> getContigLines() {
        return Collections.unmodifiableList(contigMetaData);
    }

    /**
     * check our metadata for a VCF version tag, and throw an exception if the version is out of date
     * or the version is not present
     */
    public void loadVCFVersion() {
        List<VCFHeaderLine> toRemove = new ArrayList<VCFHeaderLine>();
        for ( VCFHeaderLine line : mMetaData )
            if ( VCFHeaderVersion.isFormatString(line.getKey())) {
                toRemove.add(line);
            }
        // remove old header lines for now,
        mMetaData.removeAll(toRemove);

    }

    /**
     * load the format/info meta data maps (these are used for quick lookup by key name)
     */
    private void loadMetaDataMaps() {
        for ( VCFHeaderLine line : mMetaData ) {
            if ( line instanceof VCFInfoHeaderLine )  {
                VCFInfoHeaderLine infoLine = (VCFInfoHeaderLine)line;
                mInfoMetaData.put(infoLine.getID(), infoLine);
            } else if ( line instanceof VCFFormatHeaderLine ) {
                VCFFormatHeaderLine formatLine = (VCFFormatHeaderLine)line;
                mFormatMetaData.put(formatLine.getID(), formatLine);
            } else if ( line instanceof VCFContigHeaderLine ) {
                contigMetaData.add((VCFContigHeaderLine)line);
            } else {
                mOtherMetaData.put(line.getKey(), line);
            }
        }
    }

    /**
     * get the header fields in order they're presented in the input file (which is now required to be
     * the order presented in the spec).
     *
     * @return a set of the header fields, in order
     */
    public Set<HEADER_FIELDS> getHeaderFields() {
        return new LinkedHashSet<HEADER_FIELDS>(Arrays.asList(HEADER_FIELDS.values()));
    }

    /**
     * get the meta data, associated with this header
     *
     * @return a set of the meta data
     */
    public Set<VCFHeaderLine> getMetaData() {
        Set<VCFHeaderLine> lines = new LinkedHashSet<VCFHeaderLine>();
        lines.add(new VCFHeaderLine(VCFHeaderVersion.VCF4_0.getFormatString(), VCFHeaderVersion.VCF4_0.getVersionString()));
        lines.addAll(mMetaData);
        return Collections.unmodifiableSet(lines);
    }

    /**
     * get the genotyping sample names
     *
     * @return a list of the genotype column names, which may be empty if hasGenotypingData() returns false
     */
    public Set<String> getGenotypeSamples() {
        return mGenotypeSampleNames;
    }

    /**
     * do we have genotyping data?
     *
     * @return true if we have genotyping columns, false otherwise
     */
    public boolean hasGenotypingData() {
        return mGenotypeSampleNames.size() > 0;
    }

    /**
     * were the input samples sorted originally?
     *
     * @return true if the input samples were sorted originally, false otherwise
     */
    public boolean samplesWereAlreadySorted() {
        return samplesWereAlreadySorted;
    }

    /** @return the column count */
    public int getColumnCount() {
        return HEADER_FIELDS.values().length + (hasGenotypingData() ? mGenotypeSampleNames.size() + 1 : 0);
    }

    /**
     * @param key    the header key name
     * @return the meta data line, or null if there is none
     */
    public VCFInfoHeaderLine getInfoHeaderLine(String key) {
        return mInfoMetaData.get(key);
    }

    /**
     * @param key    the header key name
     * @return the meta data line, or null if there is none
     */
    public VCFFormatHeaderLine getFormatHeaderLine(String key) {
        return mFormatMetaData.get(key);
    }

    /**
     * @param key    the header key name
     * @return the meta data line, or null if there is none
     */
    public VCFHeaderLine getOtherHeaderLine(String key) {
        return mOtherMetaData.get(key);
    }

    /**
     * If true additional engine headers will be written to the VCF, otherwise only the walker headers will be output.
     * @return true if additional engine headers will be written to the VCF
     */
    public boolean isWriteEngineHeaders() {
        return writeEngineHeaders;
    }

    /**
     * If true additional engine headers will be written to the VCF, otherwise only the walker headers will be output.
     * @param writeEngineHeaders true if additional engine headers will be written to the VCF
     */
    public void setWriteEngineHeaders(boolean writeEngineHeaders) {
        this.writeEngineHeaders = writeEngineHeaders;
    }

    /**
     * If true, and isWriteEngineHeaders also returns true, the command line will be written to the VCF.
     * @return true if the command line will be written to the VCF
     */
    public boolean isWriteCommandLine() {
        return writeCommandLine;
    }

    /**
     * If true, and isWriteEngineHeaders also returns true, the command line will be written to the VCF.
     * @param writeCommandLine true if the command line will be written to the VCF
     */
    public void setWriteCommandLine(boolean writeCommandLine) {
        this.writeCommandLine = writeCommandLine;
    }
}
