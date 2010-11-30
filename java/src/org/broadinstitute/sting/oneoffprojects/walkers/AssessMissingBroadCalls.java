/*
 * Copyright (c) 2010.
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
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.util.*;

/**
 * Assesses calls missing from the BI set that are made by all other centers.
 * Use -B:broad,vcf and -B:1kg,vcf
 */
@Reference(window=@Window(start=-50,stop=50))
@Requires(value={})
public class AssessMissingBroadCalls extends RodWalker<Integer, Integer> {

    private static final String status_key = "BI_STATUS";
    private static final String qual_key = "BI_QUAL";

    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter writer = null;

    public void initialize() {
        final ArrayList<String> inputNames = new ArrayList<String>();
        inputNames.add("1kg");

        // setup the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit(), inputNames));
        writer.writeHeader(new VCFHeader(hInfo, new HashSet<String>()));
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return 0;

        List<Object> kgRods = tracker.getReferenceMetaData("1kg");
        // ignore places where we don't have a variant
        if ( kgRods.size() == 0 )
            return 0;

        VariantContext vc = (VariantContext)kgRods.get(0);
        Map<String, Object> attrs = new HashMap<String, Object>(vc.getAttributes());

        List<Object> biRods = tracker.getReferenceMetaData("broad");
        if ( biRods.size() == 0 )
            attrs.put(status_key, "NotCalled");
        else {
            VariantContext BIvc = (VariantContext)biRods.get(0);
            // skip the site if we called it
            if ( !BIvc.isFiltered() )
                return 0;

            attrs.put(qual_key, BIvc.getPhredScaledQual());
            
            Set<String> filters = BIvc.getFilters();
            StringBuilder sb = new StringBuilder();
            for ( String filter : filters ) {
                if ( sb.length() != 0 )
                    sb.append("-");
                sb.append(filter);
            }
            attrs.put(status_key, sb.toString());
        }

        vc = VariantContext.modifyAttributes(vc, attrs);
        writer.add(vc, ref.getBase());

        return 1;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer counter, Integer sum) {
        return counter + sum;
    }

    public void onTraversalDone(Integer sum) {}
}
