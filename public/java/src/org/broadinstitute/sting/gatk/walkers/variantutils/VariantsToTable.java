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

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.io.PrintStream;
import java.util.*;

/**
 * Emits specific fields from a VCF file to a tab-deliminated table
 *
 * <p>
 * This walker accepts a single VCF file and writes out user-selected fields from the
 * VCF as a header-containing, tab-deliminated file.  The user specifies one or more
 * fields to print with the -F NAME, each of which appears as a single column in
 * the output file, with a header named NAME, and the value of this field in the VCF
 * one per line.  NAME can be any standard VCF column (CHROM, ID, QUAL) or any binding
 * in the INFO field (AC=10).  In addition, there are specially supported values like
 * EVENTLENGTH (length of the event), TRANSITION (for SNPs), HET (count of het genotypes),
 * HOM-REF (count of homozygous reference genotypes), HOM-VAR (count of homozygous variant
 * genotypes), NO-CALL (count of no-call genotypes), TYPE (the type of event), VAR (count of
 * non-reference genotypes), NSAMPLES (number of samples), NCALLED (number of called samples),
 * GQ (from the genotype field; works only for a file with a single sample), and MULTI-ALLELIC
 * (is the record from a multi-allelic site).  Note that this tool does not support capturing any
 * GENOTYPE field values.  If a VCF record is missing a value, then the tool by
 * default throws an error, but the special value NA can be emitted instead with
 * appropriate tool arguments.
 *
 * </p>
 *
 * <h2>Input</h2>
 * <p>
 * <ul>
 *     <li>A VCF file</li>
 *     <li>A list of -F fields to write</li>
 * </ul>
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A tab-delimited file containing the values of the requested fields in the VCF file
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 *     java -jar GenomeAnalysisTK.jar \
 *     -R reference.fasta
 *     -T VariantsToTable \
 *     -V file.vcf \
 *     -F CHROM -F POS -F ID -F QUAL -F AC \
 *     -o results.table
 *
 *     would produce a file that looks like:
 *
 *     CHROM    POS ID      QUAL    AC
 *     1        10  .       50      1
 *     1        20  rs10    99      10
 *     et cetera...
 * </pre>
 *
 * @author Mark DePristo
 * @since 2010
 */
public class VariantsToTable extends RodWalker<Integer, Integer> {
    /**
     * Variants from this VCF file are used by this tool as input.
     * The file must at least contain the standard VCF header lines, but
     * can be empty (i.e., no variants are contained in the file).
     */
    @Input(fullName="variant", shortName = "V", doc="Input VCF file", required=true)
    public List<RodBinding<VariantContext>> variants;

    @Output(doc="File to which results should be written",required=true)
    protected PrintStream out;

    /**
     * -F NAME can be any standard VCF column (CHROM, ID, QUAL) or any binding in the INFO field (e.g., AC=10).
     * Note that this tool does not support capturing any GENOTYPE field values.  Note this argument
     * accepts any number of inputs.  So -F CHROM -F POS is allowed.
     */
    @Argument(fullName="fields", shortName="F", doc="The name of each field to capture for output in the table", required=true)
    public List<String> fieldsToTake = new ArrayList<String>();

    @Hidden
    @Argument(fullName="genotypeFields", shortName="GF", doc="The name of each genotype field to capture for output in the table", required=false)
    public List<String> genotypeFieldsToTake = new ArrayList<String>();
    
    /**
     * By default this tool only emits values for fields where the FILTER field is either PASS or . (unfiltered).
     * Throwing this flag will cause VariantsToTable to emit values regardless of the FILTER field value.
     */
    @Advanced
    @Argument(fullName="showFiltered", shortName="raw", doc="If provided, field values from filtered records will be included in the output", required=false)
    public boolean showFiltered = false;

    /**
     * If provided, then this tool will exit with success after this number of records have been emitted to the file.
     */
    @Advanced
    @Argument(fullName="maxRecords", shortName="M", doc="If provided, we will emit at most maxRecord records to the table", required=false)
    public int MAX_RECORDS = -1;
    int nRecords = 0;

    /**
     * By default, records with multiple ALT alleles will comprise just one line of output; note that in general this can make your resulting file
     * unreadable/malformed for certain tools like R, as the representation of multi-allelic INFO field values are often comma-separated lists
     * of values.  Using the flag will cause multi-allelic records to be split into multiple lines of output (one for each allele in the ALT field);
     * INFO field values that are not lists are copied for each of the output records while only the appropriate entry is used for lists.
     */
    @Argument(fullName="splitMultiAllelic", shortName="SMA", doc="If provided, we will split multi-allelic records into multiple lines of output", required=false)
    public boolean splitMultiAllelic = false;

    /**
     * By default, this tool throws a UserException when it encounters a field without a value in some record.  This
     * is generally useful when you mistype -F CHROM, so that you get a friendly warning about CHROM not being
     * found before the tool runs through 40M 1000G records.  However, in some cases you genuinely want to allow such
     * fields (e.g., AC not being calculated for filtered records, if included).  When provided, this argument
     * will cause VariantsToTable to write out NA values for missing fields instead of throwing an error.
     */
    @Advanced
    @Argument(fullName="allowMissingData", shortName="AMD", doc="If provided, we will not require every record to contain every field", required=false)
    public boolean ALLOW_MISSING_DATA = false;
    private final static String MISSING_DATA = "NA";

    private TreeSet<String> samples = new TreeSet<String>();

    public void initialize() {

        String genotypeHeader = "";
        if (!genotypeFieldsToTake.isEmpty()) {
            Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), variants);
            TreeSet<String> vcfSamples = new TreeSet<String>(SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE));
            samples.addAll(vcfSamples);
            
            StringBuilder sb = new StringBuilder();
            sb.append("\t");
            for (final String sample : samples) {
                for (final String gf : genotypeFieldsToTake) {
                    sb.append(sample+"."+gf+"\t");
                }
            }
            genotypeHeader = sb.toString();
        }
        // print out the header
        out.println(Utils.join("\t", fieldsToTake) + genotypeHeader);
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return 0;

        for ( VariantContext vc : tracker.getValues(variants, context.getLocation())) {
            if ( showFiltered || vc.isNotFiltered() ) {
                for ( final List<String> record : extractFields(vc, fieldsToTake, genotypeFieldsToTake, samples,
                        ALLOW_MISSING_DATA, splitMultiAllelic) )
                    out.println(Utils.join("\t", record));
            }
        }
        
        return 1;
    }

    @Override
    public boolean isDone() {
        boolean done = MAX_RECORDS != -1 && nRecords >= MAX_RECORDS;
        if ( done) logger.warn("isDone() will return true to leave after " + nRecords + " records");
        return done ;
    }

    private static final boolean isWildCard(String s) {
        return s.endsWith("*");
    }

    /**
     * Utility function that returns the list of values for each field in fields from vc.
     *
     * @param vc the VariantContext whose field values we can to capture
     * @param fields a non-null list of fields to capture from VC
     * @param genotypeFields a (possibly) null) list of fields to capture from each genotype
     * @param samples set of samples in vc, can be null in case of sites-only file
     * @param allowMissingData if false, then throws a UserException if any field isn't found in vc.  Otherwise provides a value of NA
     * @param splitMultiAllelic  if true, multiallelic variants are to be split into multiple records
     * @return List of lists of field values
     */
    private static List<List<String>> extractFields(VariantContext vc, List<String> fields, List<String> genotypeFields,
                                                    Set<String> samples, boolean allowMissingData, boolean splitMultiAllelic) {
        
        final int numRecordsToProduce = splitMultiAllelic ? vc.getAlternateAlleles().size() : 1;
        final List<List<String>> records = new ArrayList<List<String>>(numRecordsToProduce);

        int numFields = fields.size();
        final boolean addGenotypeFields = (genotypeFields != null && !genotypeFields.isEmpty() && samples != null && !samples.isEmpty());
        if (addGenotypeFields)
            numFields += genotypeFields.size()*samples.size();

        for ( int i = 0; i < numRecordsToProduce; i++ )
            records.add(new ArrayList<String>(numFields));

        for ( String field : fields ) {

            if ( splitMultiAllelic && field.equals("ALT") ) { // we need to special case the ALT field when splitting out multi-allelic records
                addFieldValue(splitAltAlleles(vc), records);
            } else if ( getters.containsKey(field) ) {
                addFieldValue(getters.get(field).get(vc), records);
            } else if ( vc.hasAttribute(field) ) {
                addFieldValue(vc.getAttribute(field, null), records);
            } else if ( isWildCard(field) ) {
                Set<String> wildVals = new HashSet<String>();
                for ( Map.Entry<String,Object> elt : vc.getAttributes().entrySet()) {
                    if ( elt.getKey().startsWith(field.substring(0, field.length() - 1)) ) {
                        wildVals.add(elt.getValue().toString());
                    }
                }

                String val = MISSING_DATA;
                if ( wildVals.size() > 0 ) {
                    List<String> toVal = new ArrayList<String>(wildVals);
                    Collections.sort(toVal);
                    val = Utils.join(",", toVal);
                }

                addFieldValue(val, records);
            } else if ( ! allowMissingData ) {
                throw new UserException(String.format("Missing field %s in vc %s at %s", field, vc.getSource(), vc));
            } else {
                addFieldValue(MISSING_DATA, records);
            }
        }

        if (addGenotypeFields) {
            for (final String sample : samples) {
                for (final String gf : genotypeFields) {
                    if (vc.hasGenotype(sample) && vc.getGenotype(sample).hasAttribute(gf))
                        addFieldValue(vc.getGenotype(sample).getAttribute(gf),records);
                    else
                        addFieldValue(MISSING_DATA, records);
                }
            }
        }
        return records;
    }

    private static void addFieldValue(Object val, List<List<String>> result) {
        final int numResultRecords = result.size();
        
        // if we're trying to create a single output record, add it
        if ( numResultRecords == 1 ) {
            result.get(0).add(val.toString());
        }
        // if this field is a list of the proper size, add the appropriate entry to each record
        else if ( (val instanceof List) && ((List)val).size() == numResultRecords ) {
            final List list = (List)val;
            for ( int i = 0; i < numResultRecords; i++ )
                result.get(i).add(list.get(i).toString());
        }
        // otherwise, add the original value to all of the records
        else {
            final String valStr = val.toString();
            for ( List<String> record : result )
                record.add(valStr);
        }
    }

    public static List<List<String>> extractFields(VariantContext vc, List<String> fields, boolean allowMissingData) {
        return extractFields(vc, fields, null, null, allowMissingData, false);
    }
    //
    // default reduce -- doesn't do anything at all
    //
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer counter, Integer sum) { return counter + sum; }
    public void onTraversalDone(Integer sum) {}

    // ----------------------------------------------------------------------------------------------------
    //
    // static system for getting values from VC by name.
    //
    // ----------------------------------------------------------------------------------------------------

    public static abstract class Getter { public abstract String get(VariantContext vc); }
    public static Map<String, Getter> getters = new HashMap<String, Getter>();

    static {
        // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
        getters.put("CHROM", new Getter() { public String get(VariantContext vc) { return vc.getChr(); } });
        getters.put("POS", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getStart()); } });
        getters.put("REF", new Getter() {
            public String get(VariantContext vc) {
                StringBuilder x = new StringBuilder();
                x.append(vc.getAlleleStringWithRefPadding(vc.getReference()));
                return x.toString();
            }
        });
        getters.put("ALT", new Getter() {
            public String get(VariantContext vc) {
                StringBuilder x = new StringBuilder();
                int n = vc.getAlternateAlleles().size();
                if ( n == 0 ) return ".";

                for ( int i = 0; i < n; i++ ) {
                    if ( i != 0 ) x.append(",");
                    x.append(vc.getAlleleStringWithRefPadding(vc.getAlternateAllele(i)));
                }
                return x.toString();
            }
        });
        getters.put("EVENTLENGTH", new Getter() { public String get(VariantContext vc) {
            int maxLength = 0;
            for ( final Allele a : vc.getAlternateAlleles() ) {
                final int length = a.length() - vc.getReference().length();
                if( Math.abs(length) > Math.abs(maxLength) ) { maxLength = length; }
            }
            return Integer.toString(maxLength);
        }});
        getters.put("QUAL", new Getter() { public String get(VariantContext vc) { return Double.toString(vc.getPhredScaledQual()); } });
        getters.put("TRANSITION", new Getter() { public String get(VariantContext vc) {
            if ( vc.isSNP() && vc.isBiallelic() )
                return VariantContextUtils.isTransition(vc) ? "1" : "0";
            else
                return "-1";
        }});
        getters.put("FILTER", new Getter() { public String get(VariantContext vc) {
            return vc.isNotFiltered() ? "PASS" : Utils.join(",", vc.getFilters()); }
        });
        getters.put("ID", new Getter() { public String get(VariantContext vc) { return vc.getID(); } });
        getters.put("HET", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getHetCount()); } });
        getters.put("HOM-REF", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getHomRefCount()); } });
        getters.put("HOM-VAR", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getHomVarCount()); } });
        getters.put("NO-CALL", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getNoCallCount()); } });
        getters.put("TYPE", new Getter() { public String get(VariantContext vc) { return vc.getType().toString(); } });
        getters.put("VAR", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getHetCount() + vc.getHomVarCount()); } });
        getters.put("NSAMPLES", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getNSamples()); } });
        getters.put("NCALLED", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getNSamples() - vc.getNoCallCount()); } });
        getters.put("MULTI-ALLELIC", new Getter() { public String get(VariantContext vc) { return Boolean.toString(vc.getAlternateAlleles().size() > 1); } });
        getters.put("GQ", new Getter() { public String get(VariantContext vc) {
            if ( vc.getNSamples() > 1 ) throw new UserException("Cannot get GQ values for multi-sample VCF");
            return String.format("%.2f", -10 * vc.getGenotype(0).getLog10PError());
        }});
    }
    
    private static Object splitAltAlleles(VariantContext vc) {
        final int numAltAlleles = vc.getAlternateAlleles().size();
        if ( numAltAlleles == 1 )
            return vc.getAlleleStringWithRefPadding(vc.getAlternateAllele(0));

        final List<String> alleles = new ArrayList<String>(numAltAlleles);
        for ( Allele allele : vc.getAlternateAlleles() )
            alleles.add(vc.getAlleleStringWithRefPadding(allele));
        return alleles;
    }
}
