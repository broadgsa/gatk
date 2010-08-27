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
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.io.PrintStream;
import java.util.*;

/**
 * Combines VCF records from different sources; supports both full merges and set unions.
 * Merge: combines multiple records into a single one; if sample names overlap then they are uniquified.
 * Union: assumes each rod represents the same set of samples (although this is not enforced); using the
 *   priority list (if provided), emits a single record instance at every position represented in the rods.
 */
//@Reference(window=@Window(start=-50,stop=50))
@Requires(value={})
public class VariantsToTable extends RodWalker<Integer, Integer> {
    @Output(doc="File to which variants should be written",required=true)
    protected PrintStream out;

    @Argument(fullName="fields", shortName="F", doc="Fields to emit from the VCF, allows any VCF field, any info field, and some meta fields like nHets", required=true)
    public String FIELDS;

    @Argument(fullName="maxRecords", shortName="M", doc="Maximum number of records to emit, if provided", required=false)
    public int MAX_RECORDS = -1;
    int nRecords = 0;

    private List<String> fieldsToTake;

    public void initialize() {
        fieldsToTake = Arrays.asList(FIELDS.toUpperCase().split(","));

        out.println(Utils.join("\t", fieldsToTake));
    }

    private static abstract class Getter { public abstract String get(VariantContext vc); }
    private static Map<String, Getter> getters = new HashMap<String, Getter>();

    static {
        // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
        getters.put("CHROM", new Getter() { public String get(VariantContext vc) { return vc.getChr(); } });
        getters.put("POS", new Getter() { public String get(VariantContext vc) { return Integer.toString(vc.getStart()); } });
        getters.put("REF", new Getter() { public String get(VariantContext vc) { return vc.getReference().toString(); } });
        getters.put("ALT", new Getter() {
            public String get(VariantContext vc) {
                StringBuilder x = new StringBuilder();
                int n = vc.getAlternateAlleles().size();
                for ( int i = 0; i < n; i++ ) {
                    if ( i != 0 ) x.append(",");
                    x.append(vc.getAlternateAllele(i).toString());
                }
                return x.toString();
            }
        });
        getters.put("QUAL", new Getter() { public String get(VariantContext vc) { return Double.toString(vc.getPhredScaledQual()); } });
        getters.put("FILTER", new Getter() { public String get(VariantContext vc) { return Utils.join(",", vc.getFilters()); } });

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
                List<String> vals = new ArrayList<String>();

                for ( String field : fieldsToTake ) {
                    String val = "UNK";

                    if ( getters.containsKey(field) ) {
                        val = getters.get(field).get(vc);
                    } else if ( vc.hasAttribute(field) ) {
                        val = vc.getAttributeAsString(field);
                    }

                    vals.add(val);
                }

                out.println(Utils.join("\t", vals));
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

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer counter, Integer sum) {
        return counter + sum;
    }

    public void onTraversalDone(Integer sum) {}
}