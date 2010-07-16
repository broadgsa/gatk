package org.broadinstitute.sting.utils.genotype.vcf;


import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.genotype.CalledGenotype;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.io.*;
import java.util.*;

/**
 * this class writes VCF files
 */
public class VCFWriter {


    // the VCF header we're storing
    private VCFHeader mHeader = null;

    // the print stream we're writting to
    private BufferedWriter mWriter;

    // were filters applied?
    private boolean filtersWereAppliedToContext = false;

    // our genotype sample fields
    private static final List<VCFGenotypeRecord> mGenotypeRecords = new ArrayList<VCFGenotypeRecord>();

    // Properties only used when using VCF4.0 encoding
    Map<String, VCFHeaderLineType> typeUsedForFormatString = new HashMap<String, VCFHeaderLineType>();
    Map<String, VCFHeaderLineType> typeUsedForInfoFields = new HashMap<String, VCFHeaderLineType>();
    Map<String, Integer> numberUsedForInfoFields = new HashMap<String, Integer>();
    Map<String, Integer> numberUsedForFormatFields = new HashMap<String, Integer>();

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

                    // Record, if line corresponds to a FORMAT field, which type will be used for writing value
                    if (line.getClass() == VCFFormatHeaderLine.class) {
                        VCFFormatHeaderLine a = (VCFFormatHeaderLine)line;
                        String key = a.getName();
                        typeUsedForFormatString.put(key,a.getType());
                        int num = a.getCount();
                        numberUsedForFormatFields.put(key,num);
                    } else if (line.getClass() == VCFInfoHeaderLine.class) {
                        VCFInfoHeaderLine a = (VCFInfoHeaderLine)line;
                        String key = a.getName();
                        typeUsedForInfoFields.put(key,a.getType());
                        int num = a.getCount();
                        numberUsedForInfoFields.put(key, num);
                    } else if (line.getClass() == VCFFilterHeaderLine.class) {
                        filtersWereAppliedToContext = true;
                    }

                mWriter.write(VCFHeader.METADATA_INDICATOR + line.toString() + "\n");
            }

            // write out the column line
            StringBuilder b = new StringBuilder();
            b.append(VCFHeader.HEADER_INDICATOR);
            for (VCFHeader.HEADER_FIELDS field : header.getHeaderFields())
                b.append(field + VCFConstants.FIELD_SEPARATOR);

            if (header.hasGenotypingData()) {
                b.append("FORMAT" + VCFConstants.FIELD_SEPARATOR);
                for (String field : header.getGenotypeSamples())
                    b.append(field + VCFConstants.FIELD_SEPARATOR);
            }
            mWriter.write(b.toString() + "\n");
            mWriter.flush();  // necessary so that writing to an output stream will work
        }
        catch (IOException e) {
            throw new RuntimeException("IOException writing the VCF header", e);
        }
    }

    /**
     * output a record to the VCF file
     *
     * @param record the record to output
     */
    public void add(VCFRecord record) {
        addRecord(record);
    }

    public void addRecord(VCFRecord record) {
        addRecord(record, VCFGenotypeWriter.VALIDATION_STRINGENCY.STRICT);
    }

    public void add(VariantContext vc, byte[] refBases) {
        if ( mHeader == null )
            throw new IllegalStateException("The VCF Header must be written before records can be added");

        String vcfString = toStringEncoding(vc, mHeader, refBases);
        try {
            mWriter.write(vcfString + "\n");
            mWriter.flush();  // necessary so that writing to an output stream will work
        } catch (IOException e) {
            throw new RuntimeException("Unable to write the VCF object to a file");
        }

    }
    /**
     * output a record to the VCF file
     *
     * @param record                the record to output
     * @param validationStringency  the validation stringency
     */
    public void addRecord(VCFRecord record, VCFGenotypeWriter.VALIDATION_STRINGENCY validationStringency) {
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

    private String toStringEncoding(VariantContext vc, VCFHeader header, byte[] refBases) {
        StringBuilder builder = new StringBuilder();

        // CHROM \t POS \t ID \t REF \t ALT \t QUAL \t FILTER \t INFO
        GenomeLoc loc = vc.getLocation();

        String contig = loc.getContig();
        long position = loc.getStart();
        String ID = vc.hasAttribute("ID") ? vc.getAttributeAsString("ID") : VCFConstants.EMPTY_ID_FIELD;


        // deal with the reference
        String referenceFromVC = new String(vc.getReference().getBases());

        double qual = vc.hasNegLog10PError() ? vc.getPhredScaledQual() : -1;
        // TODO- clean up these flags and associated code

        String filters = vc.isFiltered() ? Utils.join(";", Utils.sorted(vc.getFilters())) : (filtersWereAppliedToContext ? VCFConstants.PASSES_FILTERS_v4 : VCFConstants.UNFILTERED);

        Map<Allele, VCFGenotypeEncoding> alleleMap = new HashMap<Allele, VCFGenotypeEncoding>();
        alleleMap.put(Allele.NO_CALL, new VCFGenotypeEncoding(VCFConstants.EMPTY_ALLELE)); // convenience for lookup
        List<VCFGenotypeEncoding> vcfAltAlleles = new ArrayList<VCFGenotypeEncoding>();

        int numTrailingBases = 0, numPaddingBases = 0;
        String paddingBases = "";
        String trailingBases = "";

        // search for reference allele and find trailing and padding at the end.
         // See first if all alleles have a base encoding (ie no deletions)
        // if so, add one common base to all alleles (reference at this location)
        boolean hasBasesInAllAlleles = true;
        String refString = null;
        for ( Allele a : vc.getAlleles() ) {
            if (a.isNull() || a.isNoCall())
                hasBasesInAllAlleles = false;

            if (a.isReference())
                refString = new String(a.getBases());
        }


        // Case where there is no original allele info, e.g. when reading from VCF3.3 or when vc was produced by the GATK.
        if (!hasBasesInAllAlleles) {
            trailingBases = new String(refBases);
            numTrailingBases = 1;
            position--;
        }


        for ( Allele a : vc.getAlleles() ) {

            VCFGenotypeEncoding encoding;


            String alleleString = new String(a.getBases());
            String s = trailingBases+alleleString+paddingBases;

            encoding = new VCFGenotypeEncoding(s, true);

            // overwrite reference string by possibly padded version
            if (a.isReference())
                referenceFromVC = s;

            else {
                vcfAltAlleles.add(encoding);
            }

            alleleMap.put(a, encoding);
        }

        List<String> vcfGenotypeAttributeKeys = new ArrayList<String>();
        if ( vc.hasGenotypes() ) {
            vcfGenotypeAttributeKeys.add(VCFConstants.GENOTYPE_KEY);
            for ( String key : calcVCFGenotypeKeys(vc) ) {
                vcfGenotypeAttributeKeys.add(key);
            }
        } else if ( header.hasGenotypingData() ) {
            // this needs to be done in case all samples are no-calls
            vcfGenotypeAttributeKeys.add(VCFConstants.GENOTYPE_KEY);
        }

        String genotypeFormatString = Utils.join(VCFConstants.GENOTYPE_FIELD_SEPARATOR, vcfGenotypeAttributeKeys);

        List<VCFGenotypeRecord> genotypeObjects = new ArrayList<VCFGenotypeRecord>(vc.getGenotypes().size());
        for ( Genotype g : vc.getGenotypesSortedByName() ) {
            List<VCFGenotypeEncoding> encodings = new ArrayList<VCFGenotypeEncoding>(g.getPloidy());

            for ( Allele a : g.getAlleles() ) {
                encodings.add(alleleMap.get(a));
            }

            VCFGenotypeRecord.PHASE phasing = g.genotypesArePhased() ? VCFGenotypeRecord.PHASE.PHASED : VCFGenotypeRecord.PHASE.UNPHASED;
            VCFGenotypeRecord vcfG = new VCFGenotypeRecord(g.getSampleName(), encodings, phasing);

            for ( String key : vcfGenotypeAttributeKeys ) {
                if ( key.equals(VCFConstants.GENOTYPE_KEY) )
                    continue;


                Object val = g.hasAttribute(key) ? g.getAttribute(key) : VCFConstants.MISSING_VALUE_v4;

                // some exceptions
                if ( key.equals(VCFConstants.GENOTYPE_QUALITY_KEY) ) {
                    if ( MathUtils.compareDoubles(g.getNegLog10PError(), Genotype.NO_NEG_LOG_10PERROR) == 0 )
                        val = VCFConstants.MISSING_VALUE_v4;
                    else {
                        val = String.format(VCFConstants.DOUBLE_PRECISION_FORMAT_STRING, Math.min(g.getPhredScaledQual(), VCFConstants.MAX_GENOTYPE_QUAL));
                    }

                } else if ( key.equals(VCFConstants.DEPTH_KEY) && val == null ) {
                    ReadBackedPileup pileup = (ReadBackedPileup)g.getAttribute(CalledGenotype.READBACKEDPILEUP_ATTRIBUTE_KEY);
                    if ( pileup != null )
                        val = pileup.size();
                } else if ( key.equals(VCFConstants.GENOTYPE_FILTER_KEY) ) {
                    // VCF 4.0 key for no filters is "."
                    val = g.isFiltered() ? Utils.join(";", Utils.sorted(g.getFilters())) : VCFConstants.PASSES_FILTERS_v4;
                }


                Object newVal;
                if (typeUsedForFormatString.containsKey(key)) {
                    VCFHeaderLineType formatType = typeUsedForFormatString.get(key);
                    if (!val.getClass().equals(String.class))
                        newVal = formatType.convert(String.valueOf(val), VCFCompoundHeaderLine.SupportedHeaderLineType.FORMAT);
                    else
                        newVal = val;

                }
                else {
                    newVal = val;
                }

                if (numberUsedForFormatFields.containsKey(key)){
                    int numInFormatField = numberUsedForFormatFields.get(key);
                    if (numInFormatField>1 && val.equals(VCFConstants.MISSING_VALUE_v4)) {
                        // If we have a missing field but multiple values are expected, we need to construct new string with all fields.
                        // for example for Number =2, string has to be ".,."
                        StringBuilder v = new StringBuilder(VCFConstants.MISSING_VALUE_v4);
                        for ( int i = 1; i < numInFormatField; i++ ) {
                            v.append(",");
                            v.append(VCFConstants.MISSING_VALUE_v4);
                        }
                        newVal = v.toString();
                    }
                }
                // assume that if key is absent, given string encoding suffices.
                String outputValue = formatVCFField(key, newVal);

                if ( outputValue != null )
                    vcfG.setField(key, outputValue);
            }

            genotypeObjects.add(vcfG);
        }

        mGenotypeRecords.clear();
        mGenotypeRecords.addAll(genotypeObjects);
        // info fields
        Map<String, String> infoFields = new TreeMap<String, String>();
        for ( Map.Entry<String, Object> elt : vc.getAttributes().entrySet() ) {
            String key = elt.getKey();
            if ( key.equals("ID") )
                continue;

            String outputValue = formatVCFField(key, elt.getValue());
            if ( outputValue != null )
                infoFields.put(key, outputValue);
        }



        builder.append(contig);
        builder.append(VCFConstants.FIELD_SEPARATOR);
        builder.append(position);
        builder.append(VCFConstants.FIELD_SEPARATOR);
        builder.append(ID);
        builder.append(VCFConstants.FIELD_SEPARATOR);
        builder.append(referenceFromVC);
        builder.append(VCFConstants.FIELD_SEPARATOR);

        if ( vcfAltAlleles.size() > 0 ) {
            builder.append(vcfAltAlleles.get(0));
            for ( int i = 1; i < vcfAltAlleles.size(); i++ ) {
                builder.append(",");
                builder.append(vcfAltAlleles.get(i));
            }
        } else {
            builder.append(VCFConstants.EMPTY_ALLELE);
        }
        builder.append(VCFConstants.FIELD_SEPARATOR);


        if ( qual == -1 )
            builder.append(VCFConstants.MISSING_VALUE_v4);
        else
            builder.append(String.format(VCFConstants.DOUBLE_PRECISION_FORMAT_STRING, qual));

        builder.append(VCFConstants.FIELD_SEPARATOR);


        builder.append(filters);
        builder.append(VCFConstants.FIELD_SEPARATOR);
        builder.append(createInfoString(infoFields));

        if ( genotypeFormatString != null && genotypeFormatString.length() > 0 ) {
            addGenotypeData(builder, header, genotypeFormatString, vcfAltAlleles);
        }

        return builder.toString();


    }

    /**
     * add the genotype data
     *
     * @param builder the string builder
     * @param header  the header object
     * @param genotypeFormatString Genotype formatting string
     * @param vcfAltAlleles alternate alleles at this site
     */
    private void addGenotypeData(StringBuilder builder, VCFHeader header,
                                 String genotypeFormatString, List<VCFGenotypeEncoding>vcfAltAlleles) {
        Map<String, VCFGenotypeRecord> gMap = genotypeListToMap(mGenotypeRecords);

        StringBuffer tempStr = new StringBuffer();
        if ( header.getGenotypeSamples().size() < mGenotypeRecords.size() ) {
            for ( String sample : gMap.keySet() ) {
                if ( !header.getGenotypeSamples().contains(sample) )
                    System.err.println("Sample " + sample + " is a duplicate or is otherwise not present in the header");
                else
                    header.getGenotypeSamples().remove(sample);
            }
            throw new IllegalStateException("We have more genotype samples than the header specified; please check that samples aren't duplicated");
        }
        tempStr.append(VCFConstants.FIELD_SEPARATOR + genotypeFormatString);

        String[] genotypeFormatStrings = genotypeFormatString.split(":");

        for ( String genotype : header.getGenotypeSamples() ) {
            tempStr.append(VCFConstants.FIELD_SEPARATOR);
            if ( gMap.containsKey(genotype) ) {
                VCFGenotypeRecord rec = gMap.get(genotype);
                String genotypeString = rec.toStringEncoding(vcfAltAlleles, genotypeFormatStrings, true);

                // Override default produced genotype string when there are trailing 
                String[] genotypeStrings = genotypeString.split(":");
                int lastUsedPosition = 0;
                for (int k=genotypeStrings.length-1; k >=1; k--) {
                    // see if string represents an empty field. If not, break.
                    if (!isEmptyField(genotypeStrings[k])  ) {
                        lastUsedPosition = k;
                        break;
                    }
                }
                // now reconstruct genotypeString from 0 to lastUsedPosition
                genotypeString = Utils.join(":",genotypeStrings, 0,lastUsedPosition+1);
                tempStr.append(genotypeString);
                gMap.remove(genotype);
            } else {
                tempStr.append(VCFGenotypeRecord.stringEncodingForEmptyGenotype(genotypeFormatStrings, true));
            }
        }
        if ( gMap.size() != 0 ) {
            for ( String sample : gMap.keySet() )
                System.err.println("Sample " + sample + " is being genotyped but isn't in the header.");
            throw new IllegalStateException("We failed to use all the genotype samples; there must be an inconsistancy between the header and records");
        }

        builder.append(tempStr);
    }

    boolean isEmptyField(String field) {
        // check if given genotype field is empty, ie either ".", or ".,.", or ".,.,.", etc.
        String[] fields = field.split(",");
        boolean isEmpty = true;

        for (int k=0; k < fields.length; k++) {
            if (!fields[k].equals(VCFConstants.MISSING_VALUE_v4)) {
                isEmpty = false;
                break;
            }

        }
        return isEmpty;

    }
    /**
     * create a genotype mapping from a list and their sample names
     *
     * @param list a list of genotype samples
     * @return a mapping of the sample name to VCF genotype record
     */
    private static Map<String, VCFGenotypeRecord> genotypeListToMap(List<VCFGenotypeRecord> list) {
        Map<String, VCFGenotypeRecord> map = new HashMap<String, VCFGenotypeRecord>();
        for (int i = 0; i < list.size(); i++) {
            VCFGenotypeRecord rec = list.get(i);
            map.put(rec.getSampleName(), rec);
        }
        return map;
    }

    /**
     * create the info string
     *
     * @param infoFields a map of info fields
     * @return a string representing the infomation fields
     */
    protected String createInfoString(Map<String,String> infoFields) {
        StringBuffer info = new StringBuffer();
        boolean isFirst = true;
        for (Map.Entry<String, String> entry : infoFields.entrySet()) {
            if ( isFirst )
                isFirst = false;
            else
                info.append(VCFConstants.INFO_FIELD_SEPARATOR);

            info.append(entry.getKey());

            if ( entry.getValue() != null && !entry.getValue().equals("") ) {
                int numVals = 1;
                String key = entry.getKey();
                if (numberUsedForInfoFields.containsKey(key)) {
                    numVals = numberUsedForInfoFields.get(key);
                }

                // take care of unbounded encoding
                // TODO - workaround for "-1" in original INFO header structure
                if (numVals == VCFInfoHeaderLine.UNBOUNDED || numVals < 0)
                    numVals = 1;

                if (numVals > 0) {
                    info.append("=");
                    info.append(entry.getValue());
                }
            }
        }
        return info.length() == 0 ? VCFConstants.EMPTY_INFO_FIELD : info.toString();
    }

    private static String formatVCFField(String key, Object val) {
        String result;
        if ( val == null )
            result = VCFGenotypeRecord.getMissingFieldValue(key);
        else if ( val instanceof Double ) {
            result = String.format("%.2f", (Double)val);
        }
        else if ( val instanceof Boolean )
            result = (Boolean)val ? "" : null; // empty string for true, null for false
        else if ( val instanceof List ) {
            List list = (List)val;
            if ( list.size() == 0 )
                return formatVCFField(key, null);
            StringBuffer sb = new StringBuffer(formatVCFField(key, list.get(0)));
            for ( int i = 1; i < list.size(); i++) {
                sb.append(",");
                sb.append(formatVCFField(key, list.get(i)));
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
