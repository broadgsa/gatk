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
import org.broadinstitute.sting.utils.MathUtils;
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
 * in the INFO field (AC=10).  Note that this tool does not support capturing any
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
    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Output(doc="File to which results should be written",required=true)
    protected PrintStream out;

    /**
     * -F NAME can be any standard VCF column (CHROM, ID, QUAL) or any binding in the INFO field (e.g., AC=10).
     * Note that this tool does not support capturing any GENOTYPE field values.  Note this argument
     * accepts any number of inputs.  So -F CHROM -F POS is allowed.
     */
    @Argument(fullName="fields", shortName="F", doc="The name of each field to capture for output in the table", required=true)
    public List<String> fieldsToTake = new ArrayList<String>();

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
     * By default, only biallelic (REF=A, ALT=B) sites are including in the output.  If this flag is provided, then
     * VariantsToTable will emit field values for records with multiple ALT alleles.  Note that in general this
     * can make your resulting file unreadable and malformated according to tools like R, as the representation of
     * multi-allelic INFO field values can be lists of values.
     */
    @Advanced
     @Argument(fullName="keepMultiAllelic", shortName="KMA", doc="If provided, we will not require the site to be biallelic", required=false)
     public boolean keepMultiAllelic = false;

    @Hidden
    @Argument(fullName="logACSum", shortName="logACSum", doc="Log sum of AC instead of max value in case of multiallelic variants", required=false)
     public boolean logACSum = false;

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

    public void initialize() {
        // print out the header
        out.println(Utils.join("\t", fieldsToTake));
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return 0;

        for ( VariantContext vc : tracker.getValues(variantCollection.variants, context.getLocation())) {
            if ( (keepMultiAllelic || vc.isBiallelic()) && ( showFiltered || vc.isNotFiltered() ) ) {
                List<String> vals = extractFields(vc, fieldsToTake, ALLOW_MISSING_DATA, keepMultiAllelic, logACSum);
                out.println(Utils.join("\t", vals));
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
     * @param allowMissingData if false, then throws a UserException if any field isn't found in vc.  Otherwise
     *   provides a value of NA
     *   @param kma if true, multiallelic variants are to be kept
     *   @param logsum if true, AF and AC are computed based on sum of allele counts. Otherwise, based on allele with highest count.
     * @return
     */
    private static List<String> extractFields(VariantContext vc, List<String> fields, boolean allowMissingData, boolean kma, boolean logsum) {
        List<String> vals = new ArrayList<String>();

        for ( String field : fields ) {
            String val = "NA";

            if ( getters.containsKey(field) ) {
                val = getters.get(field).get(vc);
            } else if ( vc.hasAttribute(field) ) {
                val = vc.getAttributeAsString(field, null);
            } else if ( isWildCard(field) ) {
                Set<String> wildVals = new HashSet<String>();
                for ( Map.Entry<String,Object> elt : vc.getAttributes().entrySet()) {
                    if ( elt.getKey().startsWith(field.substring(0, field.length() - 1)) ) {
                        wildVals.add(elt.getValue().toString());
                    }
                }

                if ( wildVals.size() > 0 ) {
                    List<String> toVal = new ArrayList<String>(wildVals);
                    Collections.sort(toVal);
                    val = Utils.join(",", toVal);
                }
            } else if ( ! allowMissingData ) {
                throw new UserException(String.format("Missing field %s in vc %s at %s", field, vc.getSource(), vc));
            }

            if (field.equals("AF") || field.equals("AC")) {
                     String afo = val;

                     double af=0;
                     if (afo.contains(",")) {
                         String[] afs = afo.split(",");
                         afs[0] = afs[0].substring(1,afs[0].length());
                         afs[afs.length-1] = afs[afs.length-1].substring(0,afs[afs.length-1].length()-1);

                         double[] afd = new double[afs.length];

                         for (int k=0; k < afd.length; k++)
                             afd[k] = Double.valueOf(afs[k]);

                         if (kma && logsum)
                             af = MathUtils.sum(afd);
                         else
                         af = MathUtils.arrayMax(afd);
                         //af = Double.valueOf(afs[0]);

                     }
                     else
                         if (!afo.equals("NA"))
                             af = Double.valueOf(afo);

                val = Double.toString(af);

            }
            vals.add(val);
        }

        return vals;
    }

    public static List<String> extractFields(VariantContext vc, List<String> fields, boolean allowMissingData) {
        return extractFields(vc, fields, allowMissingData, false, false);
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
                String x = "";
                if ( vc.hasReferenceBaseForIndel() && !vc.isSNP() ) {
                    Byte refByte = vc.getReferenceBaseForIndel();
                    x=x+new String(new byte[]{refByte});
                }
                return x+vc.getReference().getDisplayString();
            }
        });
        getters.put("ALT", new Getter() {
            public String get(VariantContext vc) {
                StringBuilder x = new StringBuilder();
                int n = vc.getAlternateAlleles().size();
                if ( n == 0 ) return ".";
                if ( vc.hasReferenceBaseForIndel() && !vc.isSNP() ) {
                    Byte refByte = vc.getReferenceBaseForIndel();
                    x.append(new String(new byte[]{refByte}));
                }

                for ( int i = 0; i < n; i++ ) {
                    if ( i != 0 ) x.append(",");
                    x.append(vc.getAlternateAllele(i).getDisplayString());
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
        getters.put("ID", new Getter() { public String get(VariantContext vc) { return vc.hasID() ? vc.getID() : "."; } });
        getters.put("HET", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getHetCount()); } });
        getters.put("HOM-REF", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getHomRefCount()); } });
        getters.put("HOM-VAR", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getHomVarCount()); } });
        getters.put("NO-CALL", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getNoCallCount()); } });
        getters.put("TYPE", new Getter() { public String get(VariantContext vc) { return vc.getType().toString(); } });
        getters.put("VAR", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getHetCount() + vc.getHomVarCount()); } });
        getters.put("NSAMPLES", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getNSamples()); } });
        getters.put("NCALLED", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getNSamples() - vc.getNoCallCount()); } });
        getters.put("GQ", new Getter() { public String get(VariantContext vc) {
            if ( vc.getNSamples() > 1 ) throw new UserException("Cannot get GQ values for multi-sample VCF");
            return String.format("%.2f", 10 * vc.getGenotype(0).getNegLog10PError());
        }});
    }

}
