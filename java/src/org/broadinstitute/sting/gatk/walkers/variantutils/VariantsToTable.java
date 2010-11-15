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

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.Utils;

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
    public String FIELDS;

    @Argument(fullName="showFiltered", shortName="raw", doc="Include filtered records")
    public boolean showFiltered = false;

    @Argument(fullName="maxRecords", shortName="M", doc="Maximum number of records to emit, if provided", required=false)
    public int MAX_RECORDS = -1;
    int nRecords = 0;

    @Argument(fullName="ignoreMultiAllelic", shortName="IMA", doc="If provided, we will not require the site to be biallelic", required=false)
    public boolean ignoreMultiAllelic = false;

    private List<String> fieldsToTake;

    public void initialize() {
        fieldsToTake = Arrays.asList(FIELDS.split(","));

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
    }


    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return 0;

        if ( ++nRecords < MAX_RECORDS || MAX_RECORDS == -1 ) {
            Collection<VariantContext> vcs = tracker.getAllVariantContexts(ref, context.getLocation());
            for ( VariantContext vc : vcs) {
                if ( ! ignoreMultiAllelic || vc.isBiallelic() || ( showFiltered || vc.isNotFiltered() ) ) {
                    List<String> vals = extractFields(vc, fieldsToTake);

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

    public static List<String> extractFields(VariantContext vc, List<String> fields) {
        List<String> vals = new ArrayList<String>();

        for ( String field : fields ) {
            String val = "NA";

            if ( getters.containsKey(field) ) {
                val = getters.get(field).get(vc);
            } else if ( vc.hasAttribute(field) ) {
                val = vc.getAttributeAsString(field);
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
