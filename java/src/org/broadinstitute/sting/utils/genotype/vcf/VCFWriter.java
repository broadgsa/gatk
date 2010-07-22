package org.broadinstitute.sting.utils.genotype.vcf;


import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.CalledGenotype;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.io.*;
import java.util.*;

/**
 * this class writes VCF files
 */
public class VCFWriter {

    // the VCF header we're storing
    protected VCFHeader mHeader = null;

    // the print stream we're writting to
    protected BufferedWriter mWriter;

    // were filters applied?
    protected boolean filtersWereAppliedToContext = false;

    /**
     * create a VCF writer, given a file to write to
     *
     * @param location the file location to write to
     */
    public VCFWriter(File location) {
        FileOutputStream output;
        try {
            output = new FileOutputStream(location);
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Unable to create VCF file at location: " + location);
        }

        mWriter = new BufferedWriter(new OutputStreamWriter(output));
    }


    /**
     * create a VCF writer, given a stream to write to
     *
     * @param output   the file location to write to
     */
    public VCFWriter(OutputStream output) {
        mWriter = new BufferedWriter(new OutputStreamWriter(output));
    }

    public void writeHeader(VCFHeader header) {
        this.mHeader = header;

        try {
            // the file format field needs to be written first
            mWriter.write(VCFHeader.METADATA_INDICATOR + VCFHeaderVersion.VCF4_0.getFormatString() + "=" + VCFHeaderVersion.VCF4_0.getVersionString() + "\n");

            for ( VCFHeaderLine line : header.getMetaData() ) {
                if ( line.getKey().equals(VCFHeaderVersion.VCF4_0.getFormatString()) ||
                        line.getKey().equals(VCFHeaderVersion.VCF3_3.getFormatString()) ||
                        line.getKey().equals(VCFHeaderVersion.VCF3_2.getFormatString()) )
                    continue;

                // are the records filtered (so we know what to put in the FILTER column of passing records) ?
                if ( line instanceof VCFFilterHeaderLine )
                    filtersWereAppliedToContext = true;

                mWriter.write(VCFHeader.METADATA_INDICATOR);
                mWriter.write(line.toString());
                mWriter.write("\n");
            }

            // write out the column line
            mWriter.write(VCFHeader.HEADER_INDICATOR);
            for ( VCFHeader.HEADER_FIELDS field : header.getHeaderFields() ) {
                mWriter.write(field.toString());
                mWriter.write(VCFConstants.FIELD_SEPARATOR);
            }

            if ( header.hasGenotypingData() ) {
                mWriter.write("FORMAT");
                mWriter.write(VCFConstants.FIELD_SEPARATOR);
                for ( String sample : header.getGenotypeSamples() ) {
                   mWriter.write(sample);
                   mWriter.write(VCFConstants.FIELD_SEPARATOR);
            }                                                       }

            mWriter.write("\n");
            mWriter.flush();  // necessary so that writing to an output stream will work
        }
        catch (IOException e) {
            throw new RuntimeException("IOException writing the VCF header", e);
        }
    }

    /**
     * output a record to the VCF file
     *
     * @param record                the record to output
     */
    @Deprecated
    public void addRecord(VCFRecord record) {
        if ( mHeader == null )
            throw new IllegalStateException("The VCF Header must be written before records can be added");

        String vcfString = record.toStringEncoding(mHeader);
        try {
            mWriter.write(vcfString + "\n");
            mWriter.flush();  // necessary so that writing to an output stream will work
        } catch (IOException e) {
            throw new RuntimeException("Unable to write the VCF object to a file");
        }

    }

    /**
     * attempt to close the VCF file
     */
    public void close() {
        try {
            mWriter.flush();
            mWriter.close();
        } catch (IOException e) {
            throw new RuntimeException("Unable to close VCFFile");
        }
    }

    public void add(VariantContext vc, byte[] refBases) {
        if ( mHeader == null )
            throw new IllegalStateException("The VCF Header must be written before records can be added");
        if ( refBases == null || refBases.length < 1 )
            throw new IllegalArgumentException("The reference base must be provided to write VCF records");

        try {

            vc = VariantContextUtils.createVariantContextWithPaddedAlleles(vc, refBases);

            GenomeLoc loc = vc.getLocation();
            Map<Allele, String> alleleMap = new HashMap<Allele, String>(vc.getAlleles().size());
            alleleMap.put(Allele.NO_CALL, VCFConstants.EMPTY_ALLELE); // convenience for lookup

            // CHROM
            mWriter.write(loc.getContig());
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // POS
            mWriter.write(String.valueOf(loc.getStart()));
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // ID
            String ID = vc.hasAttribute(VariantContext.ID_KEY) ? vc.getAttributeAsString(VariantContext.ID_KEY) : VCFConstants.EMPTY_ID_FIELD;
            mWriter.write(ID);
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // REF
            alleleMap.put(vc.getReference(), "0");
            String refString = makeAlleleString(vc.getReference());
            mWriter.write(refString);
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // ALT
            if ( vc.isVariant() ) {
                Allele altAllele = vc.getAlternateAllele(0);
                alleleMap.put(altAllele, "1");
                String alt = makeAlleleString(altAllele);
                mWriter.write(alt);

                for (int i = 1; i < vc.getAlternateAlleles().size(); i++) {
                    altAllele = vc.getAlternateAllele(i);
                    alleleMap.put(altAllele, String.valueOf(i+1));
                    alt = makeAlleleString(altAllele);
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
            String filters = vc.isFiltered() ? Utils.join(";", Utils.sorted(vc.getFilters())) : (filtersWereAppliedToContext || vc.filtersWereApplied() ? VCFConstants.PASSES_FILTERS_v4 : VCFConstants.UNFILTERED);
            mWriter.write(filters);
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // INFO
            Map<String, String> infoFields = new TreeMap<String, String>();
            for ( Map.Entry<String, Object> field : vc.getAttributes().entrySet() ) {
                String key = field.getKey();
                if ( key.equals(VariantContext.ID_KEY) || key.equals(VariantContext.REFERENCE_BASE_FOR_INDEL_KEY) )
                    continue;

                String outputValue = formatVCFField(field.getValue());
                if ( outputValue != null )
                    infoFields.put(key, outputValue);
            }
            writeInfoString(infoFields);

            // FORMAT
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
                String genotypeFormatString = Utils.join(VCFConstants.GENOTYPE_FIELD_SEPARATOR, genotypeAttributeKeys);
                mWriter.write(VCFConstants.FIELD_SEPARATOR);
                mWriter.write(genotypeFormatString);

                addGenotypeData(vc, alleleMap, genotypeAttributeKeys);
            }

            mWriter.write("\n");
            mWriter.flush();  // necessary so that writing to an output stream will work
        } catch (IOException e) {
            throw new RuntimeException("Unable to write the VCF object to a file");
        }

    }

    private String getQualValue(double qual) {
        String s = String.format(VCFConstants.DOUBLE_PRECISION_FORMAT_STRING, qual);
        if ( s.endsWith(VCFConstants.DOUBLE_PRECISION_INT_SUFFIX) )
            s = s.substring(0, s.length() - VCFConstants.DOUBLE_PRECISION_INT_SUFFIX.length());
        return s;
    }

    private String makeAlleleString(Allele allele) {
        String s = new String(allele.getBases());

        return new String(allele.getBases());
    }

    /**
     * create the info string; assumes that no values are null
     *
     * @param infoFields a map of info fields
     * @throws IOException for writer
     */
    protected void writeInfoString(Map<String, String> infoFields) throws IOException {
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
                mWriter.write(g.genotypesArePhased() ? VCFConstants.PHASED : VCFConstants.UNPHASED);
                writeAllele(g.getAllele(i), alleleMap);
            }

            List<String> attrs = new ArrayList<String>(genotypeFormatKeys.size());
            for ( String key : genotypeFormatKeys ) {
                if ( key.equals(VCFConstants.GENOTYPE_KEY) )
                    continue;

                Object val = g.hasAttribute(key) ? g.getAttribute(key) : VCFConstants.MISSING_VALUE_v4;

                // some exceptions
                if ( key.equals(VCFConstants.GENOTYPE_QUALITY_KEY) ) {
                    if ( MathUtils.compareDoubles(g.getNegLog10PError(), Genotype.NO_NEG_LOG_10PERROR) == 0 )
                        val = VCFConstants.MISSING_VALUE_v4;
                    else {
                        val = getQualValue(Math.min(g.getPhredScaledQual(), VCFConstants.MAX_GENOTYPE_QUAL));
                    }
                } else if ( key.equals(VCFConstants.DEPTH_KEY) && val == null ) {
                    ReadBackedPileup pileup = (ReadBackedPileup)g.getAttribute(CalledGenotype.READBACKEDPILEUP_ATTRIBUTE_KEY);
                    if ( pileup != null )
                        val = pileup.size();
                } else if ( key.equals(VCFConstants.GENOTYPE_FILTER_KEY) ) {
                    val = g.isFiltered() ? Utils.join(";", Utils.sorted(g.getFilters())) : (g.filtersWereApplied() ? VCFConstants.PASSES_FILTERS_v4 : VCFConstants.UNFILTERED);
                }

                VCFFormatHeaderLine metaData = mHeader.getFormatHeaderLine(key);
                if ( metaData != null ) {
                    VCFHeaderLineType formatType = metaData.getType();
                    if ( !(val instanceof String) )
                        val = formatType.convert(String.valueOf(val), VCFCompoundHeaderLine.SupportedHeaderLineType.FORMAT);

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
                if ( attrs.get(i).equals(VCFConstants.MISSING_VALUE_v4) )
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

    private void writeAllele(Allele allele, Map<Allele, String> alleleMap) throws IOException {
        String encoding = alleleMap.get(allele);
        if ( encoding == null )
            throw new StingException("Allele " + allele + " is not an allele in the variant context");
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
            List list = (List)val;
            if ( list.size() == 0 )
                return formatVCFField(null);
            StringBuffer sb = new StringBuffer(formatVCFField(list.get(0)));
            for ( int i = 1; i < list.size(); i++) {
                sb.append(",");
                sb.append(formatVCFField(list.get(i)));
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
        
        return Utils.sorted(new ArrayList<String>(keys));
    }


}
