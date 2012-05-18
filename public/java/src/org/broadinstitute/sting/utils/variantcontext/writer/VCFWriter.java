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

package org.broadinstitute.sting.utils.variantcontext.writer;

import net.sf.samtools.SAMSequenceDictionary;
import org.broad.tribble.TribbleException;
import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.io.*;
import java.lang.reflect.Array;
import java.util.*;

/**
 * this class writes VCF files
 */
class VCFWriter extends IndexingVariantContextWriter {
    private final static String VERSION_LINE = VCFHeader.METADATA_INDICATOR + VCFHeaderVersion.VCF4_1.getFormatString() + "=" + VCFHeaderVersion.VCF4_1.getVersionString();

    // the print stream we're writing to
    final protected BufferedWriter mWriter;

    // should we write genotypes or just sites?
    final protected boolean doNotWriteGenotypes;

    // the VCF header we're storing
    protected VCFHeader mHeader = null;

    // were filters applied?
    protected boolean filtersWereAppliedToContext = false;

//    /**
//     * create a VCF writer, given a file to write to
//     *
//     * @param location the file location to write to
//     */
//    public StandardVCFWriter(final File location, final SAMSequenceDictionary refDict) {
//        this(location, openOutputStream(location), refDict, true, false);
//    }
//
//    public StandardVCFWriter(File location, final SAMSequenceDictionary refDict, boolean enableOnTheFlyIndexing) {
//        this(location, openOutputStream(location), refDict, enableOnTheFlyIndexing, false);
//    }
//
//    /**
//     * create a VCF writer, given a stream to write to
//     *
//     * @param output   the file location to write to
//     * @param doNotWriteGenotypes   do not write genotypes
//     */
//    public StandardVCFWriter(final OutputStream output, final SAMSequenceDictionary refDict, final boolean doNotWriteGenotypes) {
//        this(null, output, refDict, false, doNotWriteGenotypes);
//    }

    public VCFWriter(final File location, final OutputStream output, final SAMSequenceDictionary refDict, final boolean enableOnTheFlyIndexing, boolean doNotWriteGenotypes) {
        super(writerName(location, output), location, output, refDict, enableOnTheFlyIndexing);
        mWriter = new BufferedWriter(new OutputStreamWriter(getOutputStream())); // todo -- fix buffer size
        this.doNotWriteGenotypes = doNotWriteGenotypes;
    }

    // --------------------------------------------------------------------------------
    //
    // VCFWriter interface functions
    //
    // --------------------------------------------------------------------------------

    @Override
    public void writeHeader(VCFHeader header) {
        // note we need to update the mHeader object after this call because they header
        // may have genotypes trimmed out of it, if doNotWriteGenotypes is true
        mHeader = writeHeader(header, mWriter, doNotWriteGenotypes, getVersionLine(), getStreamName());

        // determine if we use filters, so we should FORCE pass the records
        // TODO -- this might not be necessary any longer as we have unfiltered, filtered, and PASS VCs
        for ( final VCFHeaderLine line : header.getMetaData() ) {
            if ( line instanceof VCFFilterHeaderLine)
                filtersWereAppliedToContext = true;
        }
    }

    public static final String getVersionLine() {
        return VERSION_LINE;
    }

    public static VCFHeader writeHeader(VCFHeader header,
                                        final Writer writer,
                                        final boolean doNotWriteGenotypes,
                                        final String versionLine,
                                        final String streamNameForError) {
        header = doNotWriteGenotypes ? new VCFHeader(header.getMetaData()) : header;
        
        try {
            // the file format field needs to be written first
            writer.write(versionLine + "\n");

            for ( VCFHeaderLine line : header.getMetaData() ) {
                if ( VCFHeaderVersion.isFormatString(line.getKey()) )
                    continue;

                writer.write(VCFHeader.METADATA_INDICATOR);
                writer.write(line.toString());
                writer.write("\n");
            }

            // write out the column line
            writer.write(VCFHeader.HEADER_INDICATOR);
            boolean isFirst = true;
            for ( VCFHeader.HEADER_FIELDS field : header.getHeaderFields() ) {
                if ( isFirst )
                    isFirst = false; // don't write out a field separator
                else
                    writer.write(VCFConstants.FIELD_SEPARATOR);
                writer.write(field.toString());
            }

            if ( header.hasGenotypingData() ) {
                writer.write(VCFConstants.FIELD_SEPARATOR);
                writer.write("FORMAT");
                for ( String sample : header.getGenotypeSamples() ) {
                    writer.write(VCFConstants.FIELD_SEPARATOR);
                    writer.write(sample);
                }
            }

            writer.write("\n");
            writer.flush();  // necessary so that writing to an output stream will work
        }
        catch (IOException e) {
            throw new ReviewedStingException("IOException writing the VCF header to " + streamNameForError, e);
        }

        return header;
    }

    /**
     * attempt to close the VCF file
     */
    @Override
    public void close() {
        // try to close the vcf stream
        try {
            mWriter.flush();
            mWriter.close();
        } catch (IOException e) {
            throw new ReviewedStingException("Unable to close " + getStreamName(), e);
        }

        super.close();
    }

    /**
     * add a record to the file
     *
     * @param vc      the Variant Context object
     */
    @Override
    public void add(VariantContext vc) {
        if ( mHeader == null )
            throw new IllegalStateException("The VCF Header must be written before records can be added: " + getStreamName());

        if ( doNotWriteGenotypes )
            vc = new VariantContextBuilder(vc).noGenotypes().make();

        try {
            vc = VariantContextUtils.createVariantContextWithPaddedAlleles(vc, false);
            super.add(vc);

            Map<Allele, String> alleleMap = buildAlleleMap(vc);

            // CHROM
            mWriter.write(vc.getChr());
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // POS
            mWriter.write(String.valueOf(vc.getStart()));
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // ID
            String ID = vc.getID();
            mWriter.write(ID);
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // REF
            String refString = vc.getReference().getDisplayString();
            mWriter.write(refString);
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // ALT
            if ( vc.isVariant() ) {
                Allele altAllele = vc.getAlternateAllele(0);
                String alt = altAllele.getDisplayString();
                mWriter.write(alt);

                for (int i = 1; i < vc.getAlternateAlleles().size(); i++) {
                    altAllele = vc.getAlternateAllele(i);
                    alt = altAllele.getDisplayString();
                    mWriter.write(",");
                    mWriter.write(alt);
                }
            } else {
                mWriter.write(VCFConstants.EMPTY_ALTERNATE_ALLELE_FIELD);
            }
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // QUAL
            if ( !vc.hasLog10PError() )
                mWriter.write(VCFConstants.MISSING_VALUE_v4);
            else
                mWriter.write(getQualValue(vc.getPhredScaledQual()));
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // FILTER
            String filters = getFilterString(vc, filtersWereAppliedToContext);
            mWriter.write(filters);
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // INFO
            Map<String, String> infoFields = new TreeMap<String, String>();
            for ( Map.Entry<String, Object> field : vc.getAttributes().entrySet() ) {
                String key = field.getKey();
                String outputValue = formatVCFField(field.getValue());
                if ( outputValue != null )
                    infoFields.put(key, outputValue);
            }
            writeInfoString(infoFields);

            // FORMAT
            final GenotypesContext gc = vc.getGenotypes();
            if ( gc instanceof LazyGenotypesContext && ((LazyGenotypesContext)gc).getUnparsedGenotypeData() != null) {
                mWriter.write(VCFConstants.FIELD_SEPARATOR);
                mWriter.write(((LazyGenotypesContext)gc).getUnparsedGenotypeData().toString());
            } else {
                List<String> genotypeAttributeKeys = new ArrayList<String>();
                if ( vc.hasGenotypes() ) {
                    genotypeAttributeKeys.addAll(calcVCFGenotypeKeys(vc));
                }

                if ( genotypeAttributeKeys.size() > 0 ) {
                    String genotypeFormatString = ParsingUtils.join(VCFConstants.GENOTYPE_FIELD_SEPARATOR, genotypeAttributeKeys);
                    mWriter.write(VCFConstants.FIELD_SEPARATOR);
                    mWriter.write(genotypeFormatString);

                    addGenotypeData(vc, alleleMap, genotypeAttributeKeys);
                }
            }
            
            mWriter.write("\n");
            mWriter.flush();  // necessary so that writing to an output stream will work
        } catch (IOException e) {
            throw new RuntimeException("Unable to write the VCF object to " + getStreamName());
        }
    }

    public static Map<Allele, String> buildAlleleMap(final VariantContext vc) {
        final Map<Allele, String> alleleMap = new HashMap<Allele, String>(vc.getAlleles().size()+1);
        alleleMap.put(Allele.NO_CALL, VCFConstants.EMPTY_ALLELE); // convenience for lookup

        final List<Allele> alleles = vc.getAlleles();
        for ( int i = 0; i < alleles.size(); i++ ) {
            alleleMap.put(alleles.get(i), String.valueOf(i));
        }

        return alleleMap;
    }

    // --------------------------------------------------------------------------------
    //
    // implementation functions
    //
    // --------------------------------------------------------------------------------

    public static final String getFilterString(final VariantContext vc) {
        return getFilterString(vc, false);
    }

    public static final String getFilterString(final VariantContext vc, boolean forcePASS) {
        return vc.isFiltered() ? ParsingUtils.join(";", ParsingUtils.sortList(vc.getFilters())) : (forcePASS || vc.filtersWereApplied() ? VCFConstants.PASSES_FILTERS_v4 : VCFConstants.UNFILTERED);
    }

    private String getQualValue(double qual) {
        String s = String.format(VCFConstants.DOUBLE_PRECISION_FORMAT_STRING, qual);
        if ( s.endsWith(VCFConstants.DOUBLE_PRECISION_INT_SUFFIX) )
            s = s.substring(0, s.length() - VCFConstants.DOUBLE_PRECISION_INT_SUFFIX.length());
        return s;
    }

    /**
     * create the info string; assumes that no values are null
     *
     * @param infoFields a map of info fields
     * @throws IOException for writer
     */
    private void writeInfoString(Map<String, String> infoFields) throws IOException {
        if ( infoFields.isEmpty() ) {
            mWriter.write(VCFConstants.EMPTY_INFO_FIELD);
            return;
        }

        boolean isFirst = true;
        for ( Map.Entry<String, String> entry : infoFields.entrySet() ) {
            if ( isFirst )
                isFirst = false;
            else
                mWriter.write(VCFConstants.INFO_FIELD_SEPARATOR);

            String key = entry.getKey();
            mWriter.write(key);

            if ( !entry.getValue().equals("") ) {
                VCFInfoHeaderLine metaData = mHeader.getInfoHeaderLine(key);
                if ( metaData == null || metaData.getCountType() != VCFHeaderLineCount.INTEGER || metaData.getCount() != 0 ) {
                    mWriter.write("=");
                    mWriter.write(entry.getValue());
                }
            }
        }
    }

    /**
     * add the genotype data
     *
     * @param vc                     the variant context
     * @param genotypeFormatKeys  Genotype formatting string
     * @param alleleMap              alleles for this context
     * @throws IOException for writer
     */
    private void addGenotypeData(VariantContext vc, Map<Allele, String> alleleMap, List<String> genotypeFormatKeys)
    throws IOException {

        for ( String sample : mHeader.getGenotypeSamples() ) {
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            Genotype g = vc.getGenotype(sample);
            if ( g == null ) {
                // TODO -- The VariantContext needs to know what the general ploidy is of the samples
                // TODO -- We shouldn't be assuming diploid genotypes here!
                mWriter.write(VCFConstants.EMPTY_GENOTYPE);
                continue;
            }

            List<String> attrs = new ArrayList<String>(genotypeFormatKeys.size());
            for ( String key : genotypeFormatKeys ) {

                if ( key.equals(VCFConstants.GENOTYPE_KEY) ) {
                    if ( !g.isAvailable() ) {
                        throw new ReviewedStingException("GTs cannot be missing for some samples if they are available for others in the record");
                    }

                    writeAllele(g.getAllele(0), alleleMap);
                    for (int i = 1; i < g.getPloidy(); i++) {
                        mWriter.write(g.isPhased() ? VCFConstants.PHASED : VCFConstants.UNPHASED);
                        writeAllele(g.getAllele(i), alleleMap);
                    }

                    continue;
                }

                Object val = g.hasAttribute(key) ? g.getAttribute(key) : VCFConstants.MISSING_VALUE_v4;

                // some exceptions
                if ( key.equals(VCFConstants.GENOTYPE_QUALITY_KEY) ) {
                    if ( ! g.hasLog10PError() )
                        val = VCFConstants.MISSING_VALUE_v4;
                    else {
                        val = getQualValue(Math.min(g.getPhredScaledQual(), VCFConstants.MAX_GENOTYPE_QUAL));
                    }
                } else if ( key.equals(VCFConstants.GENOTYPE_FILTER_KEY) ) {
                    val = g.isFiltered() ? ParsingUtils.join(";", ParsingUtils.sortList(g.getFilters())) : (g.filtersWereApplied() ? VCFConstants.PASSES_FILTERS_v4 : VCFConstants.UNFILTERED);
                }

                VCFFormatHeaderLine metaData = mHeader.getFormatHeaderLine(key);
                if ( metaData != null ) {
                    int numInFormatField = metaData.getCount(vc.getAlternateAlleles().size());
                    if ( numInFormatField > 1 && val.equals(VCFConstants.MISSING_VALUE_v4) ) {
                        // If we have a missing field but multiple values are expected, we need to construct a new string with all fields.
                        // For example, if Number=2, the string has to be ".,."
                        StringBuilder sb = new StringBuilder(VCFConstants.MISSING_VALUE_v4);
                        for ( int i = 1; i < numInFormatField; i++ ) {
                            sb.append(",");
                            sb.append(VCFConstants.MISSING_VALUE_v4);
                        }
                        val = sb.toString();
                    }
                }

                // assume that if key is absent, then the given string encoding suffices
                String outputValue = formatVCFField(val);
                if ( outputValue != null )
                    attrs.add(outputValue);
            }

            // strip off trailing missing values
            for (int i = attrs.size()-1; i >= 0; i--) {
                if ( isMissingValue(attrs.get(i)) )
                    attrs.remove(i);
                else
                    break;
            }

            for (int i = 0; i < attrs.size(); i++) {
                if ( i > 0 || genotypeFormatKeys.contains(VCFConstants.GENOTYPE_KEY) )
                    mWriter.write(VCFConstants.GENOTYPE_FIELD_SEPARATOR);
                mWriter.write(attrs.get(i));
            }
        }
    }

    private boolean isMissingValue(String s) {
        // we need to deal with the case that it's a list of missing values
        return (countOccurrences(VCFConstants.MISSING_VALUE_v4.charAt(0), s) + countOccurrences(',', s) == s.length());
    }

    private void writeAllele(Allele allele, Map<Allele, String> alleleMap) throws IOException {
        String encoding = alleleMap.get(allele);
        if ( encoding == null )
            throw new TribbleException.InternalCodecException("Allele " + allele + " is not an allele in the variant context");
        mWriter.write(encoding);
    }

    public static String formatVCFField(Object val) {
        String result;
        if ( val == null )
            result = VCFConstants.MISSING_VALUE_v4;
        else if ( val instanceof Double )
            result = String.format(VCFConstants.DOUBLE_PRECISION_FORMAT_STRING, (Double)val);
        else if ( val instanceof Boolean )
            result = (Boolean)val ? "" : null; // empty string for true, null for false
        else if ( val instanceof List ) {
            result = formatVCFField(((List)val).toArray());
        } else if ( val.getClass().isArray() ) {
            int length = Array.getLength(val);
            if ( length == 0 )
                return formatVCFField(null);
            StringBuffer sb = new StringBuffer(formatVCFField(Array.get(val, 0)));
            for ( int i = 1; i < length; i++) {
                sb.append(",");
                sb.append(formatVCFField(Array.get(val, i)));
            }
            result = sb.toString();
        } else
            result = val.toString();

        return result;
    }

    /**
     * Determine which genotype fields are in use in the genotypes in VC
     * @param vc
     * @return an ordered list of genotype fields in use in VC.  If vc has genotypes this will always include GT first
     */
    public static List<String> calcVCFGenotypeKeys(VariantContext vc) {
        Set<String> keys = new HashSet<String>();

        boolean sawGoodGT = false;
        boolean sawGoodQual = false;
        boolean sawGenotypeFilter = false;
        for ( final Genotype g : vc.getGenotypes() ) {
            keys.addAll(g.getAttributes().keySet());
            if ( g.isAvailable() )
                sawGoodGT = true;
            if ( g.hasLog10PError() )
                sawGoodQual = true;
            if (g.isFiltered() && g.isCalled())
                sawGenotypeFilter = true;
        }

        if ( sawGoodQual )
            keys.add(VCFConstants.GENOTYPE_QUALITY_KEY);

        if (sawGenotypeFilter)
            keys.add(VCFConstants.GENOTYPE_FILTER_KEY);

        List<String> sortedList = ParsingUtils.sortList(new ArrayList<String>(keys));

        // make sure the GT is first
        if ( sawGoodGT ) {
            List<String> newList = new ArrayList<String>(sortedList.size()+1);
            newList.add(VCFConstants.GENOTYPE_KEY);
            newList.addAll(sortedList);
            sortedList = newList;
        }

        if ( sortedList.isEmpty() && vc.hasGenotypes() ) {
            // this needs to be done in case all samples are no-calls
            return Collections.singletonList(VCFConstants.GENOTYPE_KEY);
        } else {
            return sortedList;
        }
    }


    private static int countOccurrences(char c, String s) {
           int count = 0;
           for (int i = 0; i < s.length(); i++) {
               count += s.charAt(i) == c ? 1 : 0;
           }
           return count;
    }
}
