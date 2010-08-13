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

package org.broadinstitute.sting.playground.gatk.walkers.vcftools;

import org.apache.commons.jexl2.Expression;
import org.apache.commons.jexl2.JexlContext;
import org.apache.commons.jexl2.MapContext;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.vcf.VCFUtils;
import org.broadinstitute.sting.utils.genotype.vcf.*;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFFilterHeaderLine;
import org.broad.tribble.vcf.VCFHeader;

import java.util.*;

/**
 * Selects variant calls for output from a user-supplied VCF file using a number of user-selectable, parameterizable criteria.  [TODO -- update to new walker style]
 */
@Requires(value={},referenceMetaData=@RMD(name="variant", type= ReferenceOrderedDatum.class))
public class VariantSelect extends RodWalker<Integer, Integer> {
    @Argument(fullName="match", shortName="match", doc="Expression used with INFO fields to select VCF records for inclusion in the output VCF(see wiki docs for more info)", required=false)
    protected String[] MATCH_STRINGS = new String[]{null};

    private VCFWriter writer = null;

    class MatchExp {
        String name;
        String expStr;
        Expression exp;

        public MatchExp(String name, String str, Expression exp) {
            this.name = name;
            this.expStr = str;
            this.exp = exp;
        }
    }

    private List<MatchExp> matchExpressions = new ArrayList<MatchExp>();

    public void initialize() {
        for ( int i = 0; i < MATCH_STRINGS.length; i++ ) {
            if ( MATCH_STRINGS[i] != null )  {
                try {
                    Expression filterExpression = VariantContextUtils.engine.createExpression(MATCH_STRINGS[i]);
                    matchExpressions.add(new MatchExp(String.format("match-%d", i), MATCH_STRINGS[i], filterExpression));
                } catch (Exception e) {
                    throw new StingException("Invalid expression used (" + MATCH_STRINGS[i] + "). Please see the JEXL docs for correct syntax.");
                }
            }
        }

        // setup the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("source", "VariantSelect"));
        hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));

        for ( MatchExp exp : matchExpressions ) {
            hInfo.add(new VCFFilterHeaderLine(exp.name, exp.expStr));
        }

        writer = new VCFWriterImpl(out);
        Set<String> samples = SampleUtils.getSampleListWithVCFHeader(getToolkit(), Arrays.asList("variant"));

        final VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        writer.writeHeader(vcfHeader);
    }

    public Integer reduceInit() { return 0; }

    /**
     * For each site of interest, rescore the genotype likelihoods by applying the specified feature set.
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return 1 if the locus was successfully processed, 0 if otherwise
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        VariantContext vc = tracker.getVariantContext(ref, "variant", null, context.getLocation(), false);
        // ignore places where we don't have a variant
        if ( vc == null )
            return 0;

        boolean someoneMatched = false;
        for ( MatchExp exp : matchExpressions ) {
            Map<String, Object> infoMap = new HashMap<String, Object>(vc.getAttributes());
            infoMap.put("QUAL", String.valueOf(vc.getPhredScaledQual()));

            JexlContext jContext = new MapContext(infoMap);

            try {
                //System.out.printf("Matching %s vs. %s%n", infoMap, exp.expStr);
                if ( (Boolean)exp.exp.evaluate(jContext) ) {
                    //System.out.printf("  => Matched%n");
                    someoneMatched = true;
                    break;
                }
            } catch (Exception e) {
                throw new StingException(e.getMessage());
            }
        }

        if ( someoneMatched )
            writer.add(vc, ref.getBase());

        return 1;
    }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {
        writer.close();
    }
}
