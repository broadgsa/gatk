/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.variant.variantcontext.VariantContextUtils;

import java.io.PrintStream;
import java.lang.reflect.Array;
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
@DocumentedGATKFeature( groupName = "Variant Evaluation and Manipulation Tools", extraDocs = {CommandLineGATK.class} )
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
     * Note that to capture GENOTYPE (FORMAT) field values, see the GF argument.  This argument accepts any number
     * of inputs.  So -F CHROM -F POS is allowed.
     */
    @Argument(fullName="fields", shortName="F", doc="The name of each field to capture for output in the table", required=false)
    public List<String> fieldsToTake = new ArrayList<String>();

    /**
     * -GF NAME can be any binding in the FORMAT field (e.g., GQ, PL).
     * Note this argument accepts any number of inputs.  So -F GQ -F PL is allowed.
     */
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
     * If provided, then this tool will exit with success after this number of VCF records have been emitted to the file.
     */
    @Argument(fullName="maxRecords", shortName="M", doc="If provided, we will emit at most maxRecord records to the table", required=false)
    public int MAX_RECORDS = -1;
    long nRecords = 0L;

    /**
     * By default, records with multiple ALT alleles will comprise just one line of output; note that in general this can make your resulting file
     * unreadable/malformed for certain tools like R, as the representation of multi-allelic INFO field values are often comma-separated lists
     * of values.  Using the flag will cause multi-allelic records to be split into multiple lines of output (one for each allele in the ALT field);
     * INFO field values that are not lists are copied for each of the output records while only the appropriate entry is used for lists.
     */
    @Argument(fullName="splitMultiAllelic", shortName="SMA", doc="If provided, we will split multi-allelic records into multiple lines of output", required=false)
    public boolean splitMultiAllelic = false;

    /**
     * By default, this tool emits one line per usable VCF record (or per allele if the -SMA flag is provided).  Using the -moltenize flag
     * will cause records to be split into multiple lines of output: one for each field provided with -F or one for each combination of sample
     * and field provided with -GF.  Note that the "Sample" column for -F fields will always be "site".
     */
    @Advanced
    @Argument(fullName="moltenize", shortName="moltenize", doc="If provided, we will produce molten output", required=false)
    public boolean moltenizeOutput = false;

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

    private final List<String> samples = new ArrayList<String>();

    public void initialize() {

        if ( !genotypeFieldsToTake.isEmpty() ) {
            Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), variants);
            TreeSet<String> vcfSamples = new TreeSet<String>(SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE));
            samples.addAll(vcfSamples);

            // optimization: if there are no samples, we don't have to worry about any genotype fields
            if ( samples.isEmpty() )
                genotypeFieldsToTake.clear();
        }

        // print out the header
        if ( moltenizeOutput ) {
            out.println("RecordID\tSample\tVariable\tValue");
        } else {
            final String baseHeader = Utils.join("\t", fieldsToTake);
            final String genotypeHeader = createGenotypeHeader(genotypeFieldsToTake, samples);
            final String separator = (!baseHeader.isEmpty() && !genotypeHeader.isEmpty()) ? "\t" : "";
            out.println(baseHeader + separator + genotypeHeader);
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return 0;

        for ( VariantContext vc : tracker.getValues(variants, context.getLocation())) {
            if ( showFiltered || vc.isNotFiltered() ) {
                nRecords++;
                for ( final List<String> record : extractFields(vc, fieldsToTake, genotypeFieldsToTake, samples, ALLOW_MISSING_DATA, splitMultiAllelic) ) {
                    if ( moltenizeOutput )
                        emitMoltenizedOutput(record);
                    else
                        out.println(Utils.join("\t", record));
                }
            }
        }
        
        return 1;
    }

    @Override
    public boolean isDone() {
        return (MAX_RECORDS != -1 && nRecords >= MAX_RECORDS);
    }

    private static final boolean isWildCard(String s) {
        return s.endsWith("*");
    }

    private static String createGenotypeHeader(final List<String> genotypeFieldsToTake, final List<String> samples) {
        boolean firstEntry = true;

        final StringBuilder sb = new StringBuilder();
        for ( final String sample : samples ) {
            for ( final String gf : genotypeFieldsToTake ) {
                if ( firstEntry )
                    firstEntry = false;
                else
                    sb.append("\t");
                // spaces in sample names are legal but wreak havoc in R data frames
                sb.append(sample.replace(" ","_"));
                sb.append(".");
                sb.append(gf);
            }
        }
        return sb.toString();
    }

    private void emitMoltenizedOutput(final List<String> record) {
        int index = 0;
        for ( final String field : fieldsToTake ) {
            out.println(String.format("%d\tsite\t%s\t%s", nRecords, field, record.get(index++)));
        }
        for ( final String sample : samples ) {
            for ( final String gf : genotypeFieldsToTake ) {
                out.println(String.format("%d\t%s\t%s\t%s", nRecords, sample.replace(" ","_"), gf, record.get(index++)));
            }
        }
    }

    /**
     * Utility function that returns the list of values for each field in fields from vc.
     *
     * @param vc                the VariantContext whose field values we can to capture
     * @param fields            a non-null list of fields to capture from VC
     * @param genotypeFields    a (possibly null) list of fields to capture from each genotype
     * @param samples           list of samples in vc
     * @param allowMissingData  if false, then throws a UserException if any field isn't found in vc.  Otherwise provides a value of NA
     * @param splitMultiAllelic if true, multiallelic variants are to be split into multiple records
     * @return List of lists of field values
     */
    private static List<List<String>> extractFields(final VariantContext vc,
                                                    final List<String> fields,
                                                    final List<String> genotypeFields,
                                                    final List<String> samples,
                                                    final boolean allowMissingData,
                                                    final boolean splitMultiAllelic) {
        
        final int numRecordsToProduce = splitMultiAllelic ? vc.getAlternateAlleles().size() : 1;
        final List<List<String>> records = new ArrayList<List<String>>(numRecordsToProduce);

        int numFields = fields.size();
        final boolean addGenotypeFields = genotypeFields != null && !genotypeFields.isEmpty();
        if ( addGenotypeFields )
            numFields += genotypeFields.size() * samples.size();

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

        if ( addGenotypeFields ) {
            for ( final String sample : samples ) {
                for ( final String gf : genotypeFields ) {
                    if ( vc.hasGenotype(sample) && vc.getGenotype(sample).hasAnyAttribute(gf) ) {
                        if ( gf.equals(VCFConstants.GENOTYPE_KEY) )
                            addFieldValue(vc.getGenotype(sample).getGenotypeString(true), records);
                        else
                            addFieldValue(vc.getGenotype(sample).getAnyAttribute(gf), records);
                    }
                    else
                        addFieldValue(MISSING_DATA, records);
                }
            }
        }

        return records;
    }

    private static void addFieldValue(final Object val, final List<List<String>> result) {
        final int numResultRecords = result.size();
        
        // if we're trying to create a single output record, add it
        if ( numResultRecords == 1 ) {
            result.get(0).add(prettyPrintObject(val));
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

    private static String prettyPrintObject(final Object val) {
        if ( val instanceof List )
            return prettyPrintObject(((List)val).toArray());

        if ( !val.getClass().isArray() )
            return val.toString();

        final int length = Array.getLength(val);
        if ( length == 0 )
            return "";

        final StringBuilder sb = new StringBuilder(prettyPrintObject(Array.get(val, 0)));
        for ( int i = 1; i < length; i++ ) {
            sb.append(",");
            sb.append(prettyPrintObject(Array.get(val, i)));
        }
        return sb.toString();
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
    public static final Map<String, Getter> getters = new HashMap<String, Getter>();

    static {
        // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
        getters.put("CHROM", new Getter() { public String get(VariantContext vc) { return vc.getChr(); } });
        getters.put("POS", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getStart()); } });
        getters.put("REF", new Getter() {
            public String get(VariantContext vc) {
                StringBuilder x = new StringBuilder();
                x.append(vc.getReference().getDisplayString());
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
                    x.append(vc.getAlternateAllele(i));
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
            return vc.getAlternateAllele(0);

        return vc.getAlternateAlleles();
    }
}
