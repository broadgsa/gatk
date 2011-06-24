/*
 * Copyright (c) 2010 The Broad Institute
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

import org.broadinstitute.sting.utils.variantcontext.MutableVariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.broadinstitute.sting.utils.IndelUtils.isInsideExtendedIndel;

/**
 * Changes annotation in the truth dataset from the filter field to the INFO field.
 *
 * @author carneiro
 * @since Mar 15, 2011
 * @help.summary Changes annotation in the truth dataset from the filter field to the INFO field.
 */

public class AnnotateTruthROD extends RodWalker<Integer, Integer> {

    @Output(doc="File to which validated variants should be written", required=true)
    protected VCFWriter vcfWriter = null;

    List<String> rodNames;

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {

        List<ReferenceOrderedDataSource> rodList = getToolkit().getRodDataSources();
        rodNames = new ArrayList<String>();

        // Initialize VCF header
        Set<VCFHeaderLine> headerLines = null;
        Set<String> samples = null;
        for (ReferenceOrderedDataSource rod : rodList) {
            rodNames.add(rod.getName());
            Map<String, VCFHeader> header = VCFUtils.getVCFHeadersFromRodPrefix(getToolkit(), rod.getName());
            samples = SampleUtils.getSampleList(header, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);
            headerLines = VCFUtils.smartMergeHeaders(header.values(), logger);
            headerLines.add(new VCFHeaderLine("source", "GenotypeAndValidate"));
        }
        if (headerLines == null|| samples == null)
            throw new UserException.BadInput("You need to provide at least one ROD file to annotate");
        vcfWriter.writeHeader(new VCFHeader(headerLines, samples));
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        int linesWritten = 0;
        // For some reason RodWalkers get map calls with null trackers
        if( tracker == null )
            return linesWritten;

        for (String rod : rodNames) {
            VariantContext vc = tracker.getVariantContext(ref, rod, null, context.getLocation(), false);
            if (!isInsideExtendedIndel(vc, ref)) {
                MutableVariantContext mvc = new MutableVariantContext(vc);
                mvc.putAttribute("GV", vc.getFilters().contains("TP") ? "T":"F");
                mvc.clearFilters();
                vcfWriter.add(mvc, ref.getBase());
                linesWritten++;
            }

        }
        return linesWritten;
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce( Integer mapValue, Integer reduceSum ) {
        return reduceSum + mapValue;
    }

    public void onTraversalDone( Integer reduceSum ) {
        logger.info(reduceSum + " lines written.");
    }
}
