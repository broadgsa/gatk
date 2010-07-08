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
 * this class writers VCF files
 */
public class VCFWriter {


    // the VCF header we're storing
    private VCFHeader mHeader = null;

    // the print stream we're writting to
    BufferedWriter mWriter;

    private boolean writingVCF40Format;
    private String PASSES_FILTERS_STRING = null;

    // our genotype sample fields
    private static final List<VCFGenotypeRecord> mGenotypeRecords = new ArrayList<VCFGenotypeRecord>();

    // Properties only used when using VCF4.0 encoding
    Map<String, VCFHeaderLineType> typeUsedForFormatString = new HashMap<String, VCFHeaderLineType>();
    Map<String, VCFHeaderLineType> typeUsedForInfoFields = new HashMap<String, VCFHeaderLineType>();
    Map<String, Integer> numberUsedForInfoFields = new HashMap<String, Integer>();
    Map<String, Integer> numberUsedForFormatFields = new HashMap<String, Integer>();

    // commonly used strings that are in the standard
    private final String FORMAT_FIELD_SEPARATOR = ":";
    private static final String GENOTYPE_FIELD_SEPARATOR = ":";
    private static final String FIELD_SEPARATOR = "\t";
    private static final String FILTER_CODE_SEPARATOR = ";";
    private static final String INFO_FIELD_SEPARATOR = ";";

    // default values
    private static final String UNFILTERED = ".";
    private static final String PASSES_FILTERS_VCF3 = "0";
    private static final String PASSES_FILTERS_VCF4 = "PASS";
    private static final String EMPTY_INFO_FIELD = ".";
    private static final String EMPTY_ID_FIELD = ".";
    private static final String EMPTY_ALLELE_FIELD = ".";
    private static final String DOUBLE_PRECISION_FORMAT_STRING = "%.2f";
    private static final String MISSING_GENOTYPE_FIELD = ".";

    /**
     * create a VCF writer, given a file to write to
     *
     * @param location the file location to write to
     */
    public VCFWriter(File location) {
        this(location, false);
    }

    public VCFWriter(File location, boolean useVCF4Format) {
        this.writingVCF40Format = useVCF4Format;
        this.PASSES_FILTERS_STRING = useVCF4Format ? PASSES_FILTERS_VCF4 : PASSES_FILTERS_VCF3;

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
        // use VCF3.3 by default
        this(output, false);
    }
    public VCFWriter(OutputStream output, boolean useVCF4Format) {
        this.writingVCF40Format = useVCF4Format;
        this.PASSES_FILTERS_STRING = useVCF4Format ? PASSES_FILTERS_VCF4 : PASSES_FILTERS_VCF3;
        mWriter = new BufferedWriter(new OutputStreamWriter(output));
    }

    public void writeHeader(VCFHeader header) {
        this.mHeader = header;

        try {
            // the file format field needs to be written first
            TreeSet<VCFHeaderLine> nonFormatMetaData = new TreeSet<VCFHeaderLine>();
            for ( VCFHeaderLine line : header.getMetaData() ) {
                if (writingVCF40Format) {
                    if ( line.getKey().equals(VCFHeaderVersion.VCF4_0.getFormatString()) ) {
                        mWriter.write(VCFHeader.METADATA_INDICATOR + VCFHeaderVersion.VCF4_0.getFormatString() + "=" + VCFHeaderVersion.VCF4_0.getVersionString() + "\n");
                    } else {
                        nonFormatMetaData.add(line);
                    }

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
                    }



                }
                else {
                    if ( line.getKey().equals(VCFHeaderVersion.VCF3_3.getFormatString()) ) {
                        mWriter.write(VCFHeader.METADATA_INDICATOR + line.toString() + "\n");
                    }
                    else if ( line.getKey().equals(VCFHeaderVersion.VCF3_2.getFormatString()) ) {
                        mWriter.write(VCFHeader.METADATA_INDICATOR + VCFHeaderVersion.VCF3_2.getFormatString() + "=" + VCFHeaderVersion.VCF3_2.getVersionString() + "\n");
                    } else {
                        nonFormatMetaData.add(line);
                    }
                }
            }

            // write the rest of the header meta-data out
            for ( VCFHeaderLine line : nonFormatMetaData )
                mWriter.write(VCFHeader.METADATA_INDICATOR + line + "\n");

            // write out the column line
            StringBuilder b = new StringBuilder();
            b.append(VCFHeader.HEADER_INDICATOR);
            for (VCFHeader.HEADER_FIELDS field : header.getHeaderFields())
                b.append(field + FIELD_SEPARATOR);

            if (header.hasGenotypingData()) {
                b.append("FORMAT" + FIELD_SEPARATOR);
                for (String field : header.getGenotypeSamples())
                    b.append(field + FIELD_SEPARATOR);
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

        if (!writingVCF40Format)
            throw new IllegalStateException("VCFWriter can only support add() method with a variant context if writing VCF4.0. Use VCFWriter(output, true) when constructing object");

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
        String ID = vc.hasAttribute("ID") ? vc.getAttributeAsString("ID") : EMPTY_ID_FIELD;


        // deal with the reference
        String referenceFromVC = new String(vc.getReference().getBases());

        double qual = vc.hasNegLog10PError() ? vc.getPhredScaledQual() : -1;
        // TODO- clean up these flags and associated code
        boolean filtersWereAppliedToContext = true;
        List<String> allowedGenotypeAttributeKeys = null;

        String filters = vc.isFiltered() ? Utils.join(";", Utils.sorted(vc.getFilters())) : (filtersWereAppliedToContext ? PASSES_FILTERS_STRING : UNFILTERED);

        Map<Allele, VCFGenotypeEncoding> alleleMap = new HashMap<Allele, VCFGenotypeEncoding>();
        alleleMap.put(Allele.NO_CALL, new VCFGenotypeEncoding(VCFGenotypeRecord.EMPTY_ALLELE)); // convenience for lookup
        List<VCFGenotypeEncoding> vcfAltAlleles = new ArrayList<VCFGenotypeEncoding>();

        int numTrailingBases = 0, numPaddingBases = 0;
        String paddingBases = "";
        String trailingBases = "";

        ArrayList<Allele> originalAlleles = (ArrayList)vc.getAttribute("ORIGINAL_ALLELE_LIST");


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

        if (originalAlleles != null) {
            // if original allele info was filled when reading from a VCF4 file,
            // determine whether there was a padding base(s) at the beginning and end.
            byte previousBase  = 0;
            int cnt=0;
            boolean firstBaseCommonInAllAlleles = true;
            Allele originalReferenceAllele = null;
            for (Allele originalAllele : originalAlleles){
                // if first base of allele is common to all of them, there may have been a common base deleted from all
                byte firstBase = originalAllele.getBases()[0];
                if (cnt > 0) {
                    if (firstBase != previousBase)
                        firstBaseCommonInAllAlleles = false;
                }
                previousBase = firstBase;
                cnt++;
                
                if (originalAllele.isReference())
                    originalReferenceAllele = originalAllele;
            }


            numTrailingBases = (firstBaseCommonInAllAlleles)? 1:0;
            position -= numTrailingBases;
            
            if (originalReferenceAllele == null)
                throw new IllegalStateException("At least one Allele must be reference");

            String originalRef = new String(originalReferenceAllele.getBases());
            numPaddingBases = originalRef.length()-refString.length()-numTrailingBases;

            if (numTrailingBases > 0) {
                trailingBases = originalRef.substring(0,numTrailingBases);
            }
            if (numPaddingBases > 0)
                paddingBases = originalRef.substring(originalRef.length()-numPaddingBases,originalRef.length());
        }
        else {
            // Case where there is no original allele info, e.g. when reading from VCF3.3 or when vc was produced by the GATK.
            if (!hasBasesInAllAlleles) {
                trailingBases = new String(refBases);
                numTrailingBases = 1;
                position--;
            }
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
            vcfGenotypeAttributeKeys.add(VCFGenotypeRecord.GENOTYPE_KEY);
            for ( String key : calcVCFGenotypeKeys(vc) ) {
                if ( allowedGenotypeAttributeKeys == null || allowedGenotypeAttributeKeys.contains(key) )
                    vcfGenotypeAttributeKeys.add(key);
            }
        }
        String genotypeFormatString = Utils.join(GENOTYPE_FIELD_SEPARATOR, vcfGenotypeAttributeKeys);

        List<VCFGenotypeRecord> genotypeObjects = new ArrayList<VCFGenotypeRecord>(vc.getGenotypes().size());
        for ( Genotype g : vc.getGenotypesSortedByName() ) {
            List<VCFGenotypeEncoding> encodings = new ArrayList<VCFGenotypeEncoding>(g.getPloidy());

            for ( Allele a : g.getAlleles() ) {
                encodings.add(alleleMap.get(a));
            }

            VCFGenotypeRecord.PHASE phasing = g.genotypesArePhased() ? VCFGenotypeRecord.PHASE.PHASED : VCFGenotypeRecord.PHASE.UNPHASED;
            VCFGenotypeRecord vcfG = new VCFGenotypeRecord(g.getSampleName(), encodings, phasing);

            for ( String key : vcfGenotypeAttributeKeys ) {
                if ( key.equals(VCFGenotypeRecord.GENOTYPE_KEY) )
                    continue;


                Object val = g.hasAttribute(key) ? g.getAttribute(key) : MISSING_GENOTYPE_FIELD;

                // some exceptions
                if ( key.equals(VCFGenotypeRecord.GENOTYPE_QUALITY_KEY) ) {
                    if ( MathUtils.compareDoubles(g.getNegLog10PError(), Genotype.NO_NEG_LOG_10PERROR) == 0 )
                        val = MISSING_GENOTYPE_FIELD;
                    else {
                        // TODO - check whether we need to saturate quality to 99 as in VCF3.3 coder. For now allow unbounded values
                        // val = Math.min(g.getPhredScaledQual(), VCFGenotypeRecord.MAX_QUAL_VALUE);
                        val = g.getPhredScaledQual();
                    }

                } else if ( key.equals(VCFGenotypeRecord.DEPTH_KEY) && val == null ) {
                    ReadBackedPileup pileup = (ReadBackedPileup)g.getAttribute(CalledGenotype.READBACKEDPILEUP_ATTRIBUTE_KEY);
                    if ( pileup != null )
                        val = pileup.size();
                } else if ( key.equals(VCFGenotypeRecord.GENOTYPE_FILTER_KEY) ) {
                    // VCF 4.0 key for no filters is "."
                    val = g.isFiltered() ? Utils.join(";", Utils.sorted(g.getFilters())) : PASSES_FILTERS_STRING;
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
                    if (numInFormatField>1 && val.equals(MISSING_GENOTYPE_FIELD)) {
                        // If we have a missing field but multiple values are expected, we need to construct new string with all fields.
                        // for example for Number =2, string has to be ".,."
                        StringBuilder v = new StringBuilder(MISSING_GENOTYPE_FIELD);
                        for ( int i = 1; i < numInFormatField; i++ ) {
                            v.append(",");
                            v.append(MISSING_GENOTYPE_FIELD);
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

            // Original alleles are not for reporting but only for internal bookkeeping
            if (key.equals("ORIGINAL_ALLELE_LIST"))
                continue;

            String outputValue = formatVCFField(key, elt.getValue());
            if ( outputValue != null )
                infoFields.put(key, outputValue);
        }



        builder.append(contig);
        builder.append(FIELD_SEPARATOR);
        builder.append(position);
        builder.append(FIELD_SEPARATOR);
        builder.append(ID);
        builder.append(FIELD_SEPARATOR);
        builder.append(referenceFromVC);
        builder.append(FIELD_SEPARATOR);

        if ( vcfAltAlleles.size() > 0 ) {
            builder.append(vcfAltAlleles.get(0));
            for ( int i = 1; i < vcfAltAlleles.size(); i++ ) {
                builder.append(",");
                builder.append(vcfAltAlleles.get(i));
            }
        } else {
            builder.append(EMPTY_ALLELE_FIELD);
        }
        builder.append(FIELD_SEPARATOR);


        if ( qual == -1 )
            builder.append(MISSING_GENOTYPE_FIELD);
        else
            builder.append(String.format(DOUBLE_PRECISION_FORMAT_STRING, qual));

        builder.append(FIELD_SEPARATOR);


        builder.append(filters);
        builder.append(FIELD_SEPARATOR);
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
        tempStr.append(FIELD_SEPARATOR + genotypeFormatString);

        String[] genotypeFormatStrings = genotypeFormatString.split(":");

        for ( String genotype : header.getGenotypeSamples() ) {
            tempStr.append(FIELD_SEPARATOR);
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
            if (!fields[k].matches(".")) {
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
                info.append(INFO_FIELD_SEPARATOR);

            info.append(entry.getKey());

            if ( entry.getValue() != null && !entry.getValue().equals("") ) {
                int numVals = 1;
                if (this.writingVCF40Format) {
                    String key = entry.getKey();
                    if (numberUsedForInfoFields.containsKey(key)) {
                        numVals = numberUsedForInfoFields.get(key);
                    }

                    // take care of unbounded encoding
                    // TODO - workaround for "-1" in original INFO header structure
                    if (numVals == VCFInfoHeaderLine.UNBOUNDED || numVals < 0)
                        numVals = 1;

                }
                if (numVals > 0) {
                    info.append("=");
                    info.append(entry.getValue());
                }
            }
        }
        return info.length() == 0 ? EMPTY_INFO_FIELD : info.toString();
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
        for ( Genotype g : vc.getGenotypes().values() ) {
            keys.addAll(g.getAttributes().keySet());
            if ( g.hasNegLog10PError() )
                sawGoodQual = true;
        }

        if ( sawGoodQual )
            keys.add(VCFGenotypeRecord.GENOTYPE_QUALITY_KEY);
        return Utils.sorted(new ArrayList<String>(keys));
    }


}
