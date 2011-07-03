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

import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.PrintStream;
import java.util.*;

/**
 * Emits specific fields as dictated by the user from one or more VCF files.
 */
@Requires(value={})
public class VariantsToTable extends RodWalker<Integer, Integer> {
    @Output(doc="File to which results should be written",required=true)
    protected PrintStream out;

    @Argument(fullName="fields", shortName="F", doc="Fields to emit from the VCF, allows any VCF field, any info field, and some meta fields like nHets", required=true)
    public ArrayList<String> fieldsToTake = new ArrayList<String>();

    @Argument(fullName="showFiltered", shortName="raw", doc="Include filtered records")
    public boolean showFiltered = false;

    @Argument(fullName="maxRecords", shortName="M", doc="Maximum number of records to emit, if provided", required=false)
    public int MAX_RECORDS = -1;
    int nRecords = 0;

    @Argument(fullName="keepMultiAllelic", shortName="KMA", doc="If provided, we will not require the site to be biallelic", required=false)
    public boolean keepMultiAllelic = false;

    @Argument(fullName="allowMissingData", shortName="AMD", doc="If provided, we will not require every record to contain every field", required=false)
    public boolean ALLOW_MISSING_DATA = false;

    public void initialize() {
        out.println(Utils.join("\t", fieldsToTake));
    }

    public static abstract class Getter { public abstract String get(VariantContext vc); }
    public static Map<String, Getter> getters = new HashMap<String, Getter>();

    static {
        // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
        getters.put("CHROM", new Getter() { public String get(VariantContext vc) { return vc.getChr(); } });
        getters.put("POS", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getStart()); } });
        getters.put("REF", new Getter() { public String get(VariantContext vc) { return vc.getReference().toString(); } });
        getters.put("ALT", new Getter() {
            public String get(VariantContext vc) {
                StringBuilder x = new StringBuilder();
                int n = vc.getAlternateAlleles().size();

                if ( n == 0 ) return ".";

                for ( int i = 0; i < n; i++ ) {
                    if ( i != 0 ) x.append(",");
                    x.append(vc.getAlternateAllele(i).toString());
                }
                return x.toString();
            }
        });
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

        getters.put("HET", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getHetCount()); } });
        getters.put("HOM-REF", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getHomRefCount()); } });
        getters.put("HOM-VAR", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getHomVarCount()); } });
        getters.put("NO-CALL", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getNoCallCount()); } });
        getters.put("VAR", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getHetCount() + vc.getHomVarCount()); } });
        getters.put("NSAMPLES", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getNSamples()); } });
        getters.put("NCALLED", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getNSamples() - vc.getNoCallCount()); } });
        getters.put("GQ", new Getter() { public String get(VariantContext vc) {
            if ( vc.getNSamples() > 1 ) throw new UserException("Cannot get GQ values for multi-sample VCF");
            return String.format("%.2f", 10 * vc.getGenotype(0).getNegLog10PError());
        }});
    }


    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return 0;

        if ( ++nRecords < MAX_RECORDS || MAX_RECORDS == -1 ) {
            Collection<VariantContext> vcs = tracker.getAllVariantContexts(ref, context.getLocation());
            for ( VariantContext vc : vcs) {
                if ( (keepMultiAllelic || vc.isBiallelic()) && ( showFiltered || vc.isNotFiltered() ) ) {
                    List<String> vals = extractFields(vc, fieldsToTake, ALLOW_MISSING_DATA);
                    out.println(Utils.join("\t", vals));
                }
            }

            return 1;
        } else {
            if ( nRecords >= MAX_RECORDS ) {
                logger.warn("Calling sys exit to leave after " + nRecords + " records");
                System.exit(0); // todo -- what's the recommend way to abort like this?
            }
            return 0;
        }
    }

    private static final boolean isWildCard(String s) {
        return s.endsWith("*");
    }

    public static List<String> extractFields(VariantContext vc, List<String> fields, boolean allowMissingData) {
        List<String> vals = new ArrayList<String>();

        for ( String field : fields ) {
            String val = "NA";

            if ( getters.containsKey(field) ) {
                val = getters.get(field).get(vc);
            } else if ( vc.hasAttribute(field) ) {
                val = vc.getAttributeAsString(field);
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

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer counter, Integer sum) {
        return counter + sum;
    }

    public void onTraversalDone(Integer sum) {}
}
