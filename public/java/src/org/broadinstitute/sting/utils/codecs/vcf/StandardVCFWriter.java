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

package org.broadinstitute.sting.utils.codecs.vcf;

import org.broad.tribble.Tribble;
import org.broad.tribble.TribbleException;
import org.broad.tribble.index.DynamicIndexCreator;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.util.LittleEndianOutputStream;
import org.broad.tribble.util.ParsingUtils;
import org.broad.tribble.util.PositionalStream;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;

import java.io.*;
import java.util.*;
import java.lang.reflect.Array;

/**
 * this class writes VCF files
 */
public class StandardVCFWriter implements VCFWriter {

    // the VCF header we're storing
    protected VCFHeader mHeader = null;

    // the print stream we're writing to
    protected BufferedWriter mWriter;
    protected PositionalStream positionalStream = null;

    // were filters applied?
    protected boolean filtersWereAppliedToContext = false;

    // should we write genotypes or just sites?
    protected boolean doNotWriteGenotypes = false;

    protected DynamicIndexCreator indexer = null;
    protected File indexFile = null;
    LittleEndianOutputStream idxStream = null;
    File location = null;

    /**         
     * create a VCF writer, given a file to write to
     *
     * @param location the file location to write to
     */
    public StandardVCFWriter(File location) {
        this(location, openOutputStream(location), true, false);
    }

    public StandardVCFWriter(File location, boolean enableOnTheFlyIndexing) {
        this(location, openOutputStream(location), enableOnTheFlyIndexing, false);
    }

    /**
     * create a VCF writer, given a stream to write to
     *
     * @param output   the file location to write to
     */
    public StandardVCFWriter(OutputStream output) {
        this(output, false);
    }

    /**
     * create a VCF writer, given a stream to write to
     *
     * @param output   the file location to write to
     * @param doNotWriteGenotypes   do not write genotypes
     */
    public StandardVCFWriter(OutputStream output, boolean doNotWriteGenotypes) {
        mWriter = new BufferedWriter(new OutputStreamWriter(output));
        this.doNotWriteGenotypes = doNotWriteGenotypes;
    }

    public StandardVCFWriter(File location, OutputStream output, boolean enableOnTheFlyIndexing, boolean doNotWriteGenotypes) {
        this.location = location;

        if ( enableOnTheFlyIndexing ) {
            indexFile = Tribble.indexFile(location);
            try {
                idxStream = new LittleEndianOutputStream(new FileOutputStream(indexFile));
                //System.out.println("Creating index on the fly for " + location);
                indexer = new DynamicIndexCreator(IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME);
                indexer.initialize(location, indexer.defaultBinSize());
                positionalStream = new PositionalStream(output);
                output = positionalStream;
            } catch ( IOException ex ) {
                // No matter what we keep going, since we don't care if we can't create the index file
            }
        }

        //mWriter = new BufferedWriter(new OutputStreamWriter(new PositionalStream(output)));
        mWriter = new BufferedWriter(new OutputStreamWriter(output));
        this.doNotWriteGenotypes = doNotWriteGenotypes;
    }

    public void writeHeader(VCFHeader header) {
        mHeader = doNotWriteGenotypes ? new VCFHeader(header.getMetaData()) : header;
        
        try {
            // the file format field needs to be written first
            mWriter.write(VCFHeader.METADATA_INDICATOR + VCFHeaderVersion.VCF4_0.getFormatString() + "=" + VCFHeaderVersion.VCF4_0.getVersionString() + "\n");

            for ( VCFHeaderLine line : mHeader.getMetaData() ) {
                if ( line.getKey().equals(VCFHeaderVersion.VCF4_0.getFormatString()) ||
                        line.getKey().equals(VCFHeaderVersion.VCF3_3.getFormatString()) ||
                        line.getKey().equals(VCFHeaderVersion.VCF3_2.getFormatString()) )
                    continue;

                // are the records filtered (so we know what to put in the FILTER column of passing records) ?
                if ( line instanceof VCFFilterHeaderLine)
                    filtersWereAppliedToContext = true;

                mWriter.write(VCFHeader.METADATA_INDICATOR);
                mWriter.write(line.toString());
                mWriter.write("\n");
            }

            // write out the column line
            mWriter.write(VCFHeader.HEADER_INDICATOR);
            for ( VCFHeader.HEADER_FIELDS field : mHeader.getHeaderFields() ) {
                mWriter.write(field.toString());
                mWriter.write(VCFConstants.FIELD_SEPARATOR);
            }

            if ( mHeader.hasGenotypingData() ) {
                mWriter.write("FORMAT");
                for ( String sample : mHeader.getGenotypeSamples() ) {
                    mWriter.write(VCFConstants.FIELD_SEPARATOR);
                    mWriter.write(sample);
                }
            }

            mWriter.write("\n");
            mWriter.flush();  // necessary so that writing to an output stream will work
        }
        catch (IOException e) {
            throw new TribbleException("IOException writing the VCF header to " + locationString(), e);
        }
    }

    private String locationString() {
        return location == null ? mWriter.toString() : location.getAbsolutePath();
    }

    /**
     * attempt to close the VCF file
     */
    public void close() {
        // try to close the vcf stream
        try {
            mWriter.flush();
            mWriter.close();
        } catch (IOException e) {
            throw new TribbleException("Unable to close " + locationString() + " because of " + e.getMessage());
        }

        // try to close the index stream (keep it separate to help debugging efforts)
        if ( indexer != null ) {
            try {
                Index index = indexer.finalizeIndex(positionalStream.getPosition());
                index.write(idxStream);
                idxStream.close();
            } catch (IOException e) {
                throw new TribbleException("Unable to close index for " + locationString() + " because of " + e.getMessage());
            }
        }
    }

    protected static OutputStream openOutputStream(File location) {
        try {
            return new FileOutputStream(location);
        } catch (FileNotFoundException e) {
            throw new TribbleException("Unable to create VCF file at location: " + location);
        }
    }

    /**
     * add a record to the file
     *
     * @param vc      the Variant Context object
     * @param refBase the ref base used for indels
     */
    public void add(VariantContext vc, byte refBase) {
        add(vc, refBase, false);
    }

    /**
     * add a record to the file
     *
     * @param vc      the Variant Context object
     * @param refBase the ref base used for indels
     * @param refBaseShouldBeAppliedToEndOfAlleles *** THIS SHOULD BE FALSE EXCEPT FOR AN INDEL AT THE EXTREME BEGINNING OF A CONTIG (WHERE THERE IS NO PREVIOUS BASE, SO WE USE THE BASE AFTER THE EVENT INSTEAD)
     */
    public void add(VariantContext vc, byte refBase, boolean refBaseShouldBeAppliedToEndOfAlleles) {
        if ( mHeader == null )
            throw new IllegalStateException("The VCF Header must be written before records can be added: " + locationString());

        if ( doNotWriteGenotypes )
            vc = VariantContext.modifyGenotypes(vc, null);

        try {
            vc = VariantContext.createVariantContextWithPaddedAlleles(vc, refBase, refBaseShouldBeAppliedToEndOfAlleles);

            // if we are doing on the fly indexing, add the record ***before*** we write any bytes 
            if ( indexer != null ) indexer.addFeature(vc, positionalStream.getPosition());

            Map<Allele, String> alleleMap = new HashMap<Allele, String>(vc.getAlleles().size());
            alleleMap.put(Allele.NO_CALL, VCFConstants.EMPTY_ALLELE); // convenience for lookup

            // CHROM
            mWriter.write(vc.getChr());
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // POS
            mWriter.write(String.valueOf(vc.getStart()));
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // ID
            String ID = vc.hasID() ? vc.getID() : VCFConstants.EMPTY_ID_FIELD;
            mWriter.write(ID);
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // REF
            alleleMap.put(vc.getReference(), "0");
            String refString = vc.getReference().getDisplayString();
            mWriter.write(refString);
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // ALT
            if ( vc.isVariant() ) {
                Allele altAllele = vc.getAlternateAllele(0);
                alleleMap.put(altAllele, "1");
                String alt = altAllele.getDisplayString();
                mWriter.write(alt);

                for (int i = 1; i < vc.getAlternateAlleles().size(); i++) {
                    altAllele = vc.getAlternateAllele(i);
                    alleleMap.put(altAllele, String.valueOf(i+1));
                    alt = altAllele.getDisplayString();
                    mWriter.write(",");
                    mWriter.write(alt);
                }
            } else {
                mWriter.write(VCFConstants.EMPTY_ALTERNATE_ALLELE_FIELD);
            }
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // QUAL
            if ( !vc.hasNegLog10PError() )
                mWriter.write(VCFConstants.MISSING_VALUE_v4);
            else
                mWriter.write(getQualValue(vc.getPhredScaledQual()));
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // FILTER
            String filters = vc.isFiltered() ? ParsingUtils.join(";", ParsingUtils.sortList(vc.getFilters())) : (filtersWereAppliedToContext || vc.filtersWereApplied() ? VCFConstants.PASSES_FILTERS_v4 : VCFConstants.UNFILTERED);
            mWriter.write(filters);
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // INFO
            Map<String, String> infoFields = new TreeMap<String, String>();
            for ( Map.Entry<String, Object> field : vc.getAttributes().entrySet() ) {
                String key = field.getKey();
                if ( key.equals(VariantContext.ID_KEY) || key.equals(VariantContext.REFERENCE_BASE_FOR_INDEL_KEY) || key.equals(VariantContext.UNPARSED_GENOTYPE_MAP_KEY) || key.equals(VariantContext.UNPARSED_GENOTYPE_PARSER_KEY) )
                    continue;

                String outputValue = formatVCFField(field.getValue());
                if ( outputValue != null )
                    infoFields.put(key, outputValue);
            }
            writeInfoString(infoFields);

            // FORMAT
            if ( vc.hasAttribute(VariantContext.UNPARSED_GENOTYPE_MAP_KEY) ) {
                mWriter.write(VCFConstants.FIELD_SEPARATOR);
                mWriter.write(vc.getAttributeAsString(VariantContext.UNPARSED_GENOTYPE_MAP_KEY, ""));
            } else {
                List<String> genotypeAttributeKeys = new ArrayList<String>();
                if ( vc.hasGenotypes() ) {
                    genotypeAttributeKeys.add(VCFConstants.GENOTYPE_KEY);
                    for ( String key : calcVCFGenotypeKeys(vc) ) {
                        genotypeAttributeKeys.add(key);
                    }
                } else if ( mHeader.hasGenotypingData() ) {
                    // this needs to be done in case all samples are no-calls
                    genotypeAttributeKeys.add(VCFConstants.GENOTYPE_KEY);
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
            throw new RuntimeException("Unable to write the VCF object to " + locationString());
        }

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
                int numVals = 1;
                VCFInfoHeaderLine metaData = mHeader.getInfoHeaderLine(key);
                if ( metaData != null )
                    numVals = metaData.getCount();

                // take care of unbounded encoding
                if ( numVals == VCFInfoHeaderLine.UNBOUNDED )
                    numVals = 1;

                if ( numVals > 0 ) {
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

            writeAllele(g.getAllele(0), alleleMap);
            for (int i = 1; i < g.getPloidy(); i++) {
                mWriter.write(g.isPhased() ? VCFConstants.PHASED : VCFConstants.UNPHASED);
                writeAllele(g.getAllele(i), alleleMap);
            }

            List<String> attrs = new ArrayList<String>(genotypeFormatKeys.size());
            for ( String key : genotypeFormatKeys ) {
                if ( key.equals(VCFConstants.GENOTYPE_KEY) )
                    continue;

                Object val = g.hasAttribute(key) ? g.getAttribute(key) : VCFConstants.MISSING_VALUE_v4;

                // some exceptions
                if ( key.equals(VCFConstants.GENOTYPE_QUALITY_KEY) ) {
                    if ( Math.abs(g.getNegLog10PError() - Genotype.NO_NEG_LOG_10PERROR) < 1e-6)
                        val = VCFConstants.MISSING_VALUE_v4;
                    else {
                        val = getQualValue(Math.min(g.getPhredScaledQual(), VCFConstants.MAX_GENOTYPE_QUAL));
                    }
                } else if ( key.equals(VCFConstants.GENOTYPE_FILTER_KEY) ) {
                    val = g.isFiltered() ? ParsingUtils.join(";", ParsingUtils.sortList(g.getFilters())) : (g.filtersWereApplied() ? VCFConstants.PASSES_FILTERS_v4 : VCFConstants.UNFILTERED);
                }

                VCFFormatHeaderLine metaData = mHeader.getFormatHeaderLine(key);
                if ( metaData != null ) {
                    int numInFormatField = metaData.getCount();
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

            for (String s : attrs ) {
                mWriter.write(VCFConstants.GENOTYPE_FIELD_SEPARATOR);
                mWriter.write(s);
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

    private static String formatVCFField(Object val) {
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

    private static List<String> calcVCFGenotypeKeys(VariantContext vc) {
        Set<String> keys = new HashSet<String>();

        boolean sawGoodQual = false;
        boolean sawGenotypeFilter = false;
        for ( Genotype g : vc.getGenotypes().values() ) {
            keys.addAll(g.getAttributes().keySet());
            if ( g.hasNegLog10PError() )
                sawGoodQual = true;
            if (g.isFiltered() && g.isCalled())
                sawGenotypeFilter = true;
        }

        if ( sawGoodQual )
            keys.add(VCFConstants.GENOTYPE_QUALITY_KEY);

        if (sawGenotypeFilter)
            keys.add(VCFConstants.GENOTYPE_FILTER_KEY);

        return ParsingUtils.sortList(new ArrayList<String>(keys));
    }


    public static int countOccurrences(char c, String s) {
           int count = 0;
           for (int i = 0; i < s.length(); i++) {
               count += s.charAt(i) == c ? 1 : 0;
           }
           return count;
    }

}
