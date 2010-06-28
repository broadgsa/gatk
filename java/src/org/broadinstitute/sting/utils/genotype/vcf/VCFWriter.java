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

    // our genotype sample fields
    private static final List<VCFGenotypeRecord> mGenotypeRecords = new ArrayList<VCFGenotypeRecord>();

    // commonly used strings that are in the standard
    private final String FORMAT_FIELD_SEPARATOR = ":";
    private static final String GENOTYPE_FIELD_SEPARATOR = ":";
    private static final String FIELD_SEPARATOR = "\t";
    private static final String FILTER_CODE_SEPARATOR = ";";
    private static final String INFO_FIELD_SEPARATOR = ";";

    // default values
    private static final String UNFILTERED = ".";
    private static final String PASSES_FILTERS = "0";
    private static final String PASSES_FILTERS_VCF_4_0 = "PASS";
    private static final String EMPTY_INFO_FIELD = ".";
    private static final String EMPTY_ID_FIELD = ".";
    private static final String EMPTY_ALLELE_FIELD = ".";
    private static final String DOUBLE_PRECISION_FORMAT_STRING = "%.2f";
    private static final String MISSING_GENOTYPE_QUALITY = ".";

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
        // use 3.3 format by default
        // TODO - update to true once 4.0 is on by default
        writeHeader(header, false);
    }
    public void writeHeader(VCFHeader header, boolean useVCF4Format) {
        this.mHeader = header;
        try {
            // the file format field needs to be written first
            TreeSet<VCFHeaderLine> nonFormatMetaData = new TreeSet<VCFHeaderLine>();
            for ( VCFHeaderLine line : header.getMetaData() ) {
                if (useVCF4Format) {
                    if ( line.getKey().equals(VCFHeaderVersion.VCF4_0.getFormatString()) ) {
                        mWriter.write(VCFHeader.METADATA_INDICATOR + VCFHeaderVersion.VCF4_0.getFormatString() + "=" + VCFHeaderVersion.VCF4_0.getVersionString() + "\n");
                    } else {
                        nonFormatMetaData.add(line);
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
            for (VCFHeader.HEADER_FIELDS field : header.getHeaderFields()) b.append(field + FIELD_SEPARATOR);
            if (header.hasGenotypingData()) {
                b.append("FORMAT" + FIELD_SEPARATOR);
                for (String field : header.getGenotypeSamples()) b.append(field + FIELD_SEPARATOR);
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
        //String referenceBases = new String(vc.getReference().getBases());

        String contig = loc.getContig();
        long position = loc.getStart();
        String ID = vc.hasAttribute("ID") ? vc.getAttributeAsString("ID") : EMPTY_ID_FIELD;

        builder.append(contig);
        builder.append(FIELD_SEPARATOR);
        builder.append(position);
        builder.append(FIELD_SEPARATOR);
        builder.append(ID);
        builder.append(FIELD_SEPARATOR);

        // deal with the reference
        String referenceBases = new String(vc.getReference().getBases());

        double qual = vc.hasNegLog10PError() ? vc.getPhredScaledQual() : -1;
        // TODO- clean up these flags and associated code
        boolean filtersWereAppliedToContext = true;
        List<String> allowedGenotypeAttributeKeys = null;
        boolean filtersWereAppliedToGenotypes = false;
        String filters = vc.isFiltered() ? Utils.join(";", Utils.sorted(vc.getFilters())) : (filtersWereAppliedToContext ? PASSES_FILTERS_VCF_4_0 : UNFILTERED);

        String referenceString = new String(vc.getReference().getBases());
        Map<Allele, VCFGenotypeEncoding> alleleMap = new HashMap<Allele, VCFGenotypeEncoding>();
        alleleMap.put(Allele.NO_CALL, new VCFGenotypeEncoding(VCFGenotypeRecord.EMPTY_ALLELE)); // convenience for lookup
        List<VCFGenotypeEncoding> vcfAltAlleles = new ArrayList<VCFGenotypeEncoding>();

        boolean isComplexEvent = false;
        if (referenceBases.length()>1) {
            // complex event: by VCF 4.0, reference from previous position is included.
            position--;
            isComplexEvent = true;
        }
        for ( Allele a : vc.getAlleles() ) {

            VCFGenotypeEncoding encoding;


            String alleleString = new String(a.getBases());
            if (isComplexEvent) {// vc.getType() == VariantContext.Type.MIXED ) {
                // Mixed variants (e.g. microsatellites) have the reference before the event included
                String s = new String(refBases)+alleleString;
                encoding = new VCFGenotypeEncoding(s, true);
                if (a.isReference())
                    referenceString = s;


            } else if ( vc.getType() == VariantContext.Type.INDEL ) {
                if ( a.isNull() ) {
                    if ( a.isReference() ) {
                        // ref, where alt is insertion
                        encoding = new VCFGenotypeEncoding(new String(refBases));
                    } else {
                        // non-ref deletion
                        encoding = new VCFGenotypeEncoding(".");
                    }
                } else if ( a.isReference() ) {
                    // ref, where alt is deletion
                    encoding = new VCFGenotypeEncoding(new String(refBases));
                } else {
                    // non-ref insertion
                    referenceString = new String(refBases)+alleleString;
                    encoding = new VCFGenotypeEncoding(referenceString, true);
                }
            } else {
                // no variation, ref or alt for snp
                encoding = new VCFGenotypeEncoding(alleleString);
            }

            if ( a.isNonReference() ) {
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
            if ( filtersWereAppliedToGenotypes )
                vcfGenotypeAttributeKeys.add(VCFGenotypeRecord.GENOTYPE_FILTER_KEY);
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

                Object val = g.getAttribute(key);
                // some exceptions
                if ( key.equals(VCFGenotypeRecord.GENOTYPE_QUALITY_KEY) ) {
                    if ( MathUtils.compareDoubles(g.getNegLog10PError(), Genotype.NO_NEG_LOG_10PERROR) == 0 )
                        val = VCFGenotypeRecord.MISSING_GENOTYPE_QUALITY;
                    else
                        val = Math.min(g.getPhredScaledQual(), VCFGenotypeRecord.MAX_QUAL_VALUE);
                } else if ( key.equals(VCFGenotypeRecord.DEPTH_KEY) && val == null ) {
                    ReadBackedPileup pileup = (ReadBackedPileup)g.getAttribute(CalledGenotype.READBACKEDPILEUP_ATTRIBUTE_KEY);
                    if ( pileup != null )
                        val = pileup.size();
                } else if ( key.equals(VCFGenotypeRecord.GENOTYPE_FILTER_KEY) ) {
                    val = g.isFiltered() ? Utils.join(";", Utils.sorted(g.getFilters())) : PASSES_FILTERS;
                }

                String outputValue = formatVCFField(key, val);
                if ( outputValue != null )
                    vcfG.setField(key, outputValue);
            }

            genotypeObjects.add(vcfG);
        }

        mGenotypeRecords.clear();
        mGenotypeRecords.addAll(genotypeObjects);
        // info fields
        Map<String, String> infoFields = new HashMap<String, String>();
        for ( Map.Entry<String, Object> elt : vc.getAttributes().entrySet() ) {
            String key = elt.getKey();
            if ( key.equals("ID") )
                continue;

            String outputValue = formatVCFField(key, elt.getValue());
            if ( outputValue != null )
                infoFields.put(key, outputValue);
        }




        builder.append(referenceString);
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
            builder.append(MISSING_GENOTYPE_QUALITY);
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
     */
    private static void addGenotypeData(StringBuilder builder, VCFHeader header,
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
                tempStr.append(rec.toStringEncoding(vcfAltAlleles, genotypeFormatStrings));
                gMap.remove(genotype);
            } else {
                tempStr.append(VCFGenotypeRecord.stringEncodingForEmptyGenotype(genotypeFormatStrings));
            }
        }
        if ( gMap.size() != 0 ) {
            for ( String sample : gMap.keySet() )
                System.err.println("Sample " + sample + " is being genotyped but isn't in the header.");
            throw new IllegalStateException("We failed to use all the genotype samples; there must be an inconsistancy between the header and records");
        }

        builder.append(tempStr);
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
     * @return a string representing the infomation fields
     */
    static protected String createInfoString(Map<String,String> infoFields) {
        StringBuffer info = new StringBuffer();
        boolean isFirst = true;
        for (Map.Entry<String, String> entry : infoFields.entrySet()) {
            if ( isFirst )
                isFirst = false;
            else
                info.append(INFO_FIELD_SEPARATOR);
            info.append(entry.getKey());
            if ( entry.getValue() != null && !entry.getValue().equals("") ) {
                info.append("=");
                info.append(entry.getValue());
            }
        }
        return info.length() == 0 ? EMPTY_INFO_FIELD : info.toString();
    }

    private static String formatVCFField(String key, Object val) {
        String result;
        if ( val == null )
            result = VCFGenotypeRecord.getMissingFieldValue(key);
        else if ( val instanceof Double )
            result = String.format("%.2f", (Double)val);
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
