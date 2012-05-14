/*
 * Copyright (c) 2012, The Broad Institute
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

import org.apache.commons.io.FilenameUtils;
import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.sting.utils.text.ListFileUtils;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.io.File;
import java.util.*;

/**
 * Selects headers from a VCF source.
 * <p/>
 * <p>
 * Often, a VCF containing many headers will need to be subset in order to facilitate certain formatting guidelines.
 * SelectHeaders can be used for this purpose. Given a single VCF file, one or more headers can be extracted from the
 * file (based on a complete header name or a pattern match).
 * <p/>
 * <h2>Input</h2>
 * <p>
 * A set of VCFs.
 * </p>
 * <p/>
 * <h2>Output</h2>
 * <p>
 * A header selected VCF.
 * </p>
 * <p/>
 * <h2>Examples</h2>
 * <pre>
 * Select only the FILTER, FORMAT, and INFO headers:
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectHeaders \
 *   --variant input.vcf \
 *   -o output.vcf \
 *   -hn FILTER \
 *   -hn FORMAT \
 *   -hn INFO
 *
 * Select only the FILTER, FORMAT, and INFO headers and add in the reference file names:
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectHeaders \
 *   --variant input.vcf \
 *   -o output.vcf \
 *   -hn FILTER \
 *   -hn FORMAT \
 *   -hn INFO \
 *   -irn \
 *   -iln
 *
 * Select only the FILTER, FORMAT, and INFO headers, plus any headers with SnpEff:
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectHeaders \
 *   --variant input.vcf \
 *   -o output.vcf \
 *   -hn FILTER \
 *   -hn FORMAT \
 *   -hn INFO \
 *   -he '.*SnpEff.*'
 * </pre>
 */
@SuppressWarnings("unused")
public class SelectHeaders extends RodWalker<Integer, Integer> implements TreeReducible<Integer> {
    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Output(doc = "File to which variants should be written", required = true)
    protected VariantContextWriter vcfWriter;

    @Argument(fullName = "header_name", shortName = "hn", doc = "Include header. Can be specified multiple times", required = false)
    public Set<String> headerNames;

    @Argument(fullName = "header_expression", shortName = "he", doc = "Regular expression to select many headers from the tracks provided. Can be specified multiple times", required = false)
    public Set<String> headerExpressions;

    /**
     * Note that header exclusion takes precedence over inclusion, so that if a header is in both lists it will be excluded.
     */
    @Argument(fullName = "exclude_header_name", shortName = "xl_hn", doc = "Exclude header. Can be specified multiple times", required = false)
    public Set<String> XLheaderNames;

    /**
     * Note that reference inclusion takes precedence over other header matching. If set other reference lines may be excluded but the file name will still be added.
     */
    @Argument(fullName = "include_reference_name", shortName = "irn", doc = "If set the reference file name minus the file extension will be added to the headers", required = false)
    public boolean includeReference;

    /**
     * Note that interval name inclusion takes precedence over other header matching. If set other interval lines may be excluded but the intervals will still be added.
     */
    @Argument(fullName = "include_interval_names", shortName = "iln", doc = "If set the interval file name minus the file extension, or the command line intervals, will be added to the headers", required = false)
    public boolean includeIntervals;

    /**
     * Note that engine header inclusion takes precedence over other header matching. If set other engine lines may be excluded but the intervals will still be added.
     */
    @Hidden // TODO: Determine if others find this valuable and either remove @Hidden or remove -ieh.
    @Argument(fullName = "include_engine_headers", shortName = "ieh", doc = "If set the headers normally output by the engine will be added to the headers", required = false)
    public boolean includeEngineHeaders;

    private static final ListFileUtils.StringConverter<VCFHeaderLine> headerKey = new ListFileUtils.StringConverter<VCFHeaderLine>() {
        @Override
        public String convert(VCFHeaderLine value) {
            return value.getKey();
        }
    };

    /**
     * Set up the VCF writer, the header expressions and regexps
     */
    @Override
    public void initialize() {
        // Get list of samples to include in the output
        List<String> rodNames = Arrays.asList(variantCollection.variants.getName());

        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), logger);

        headerLines.add(new VCFHeaderLine(VCFHeader.SOURCE_KEY, "SelectHeaders"));

        // Select only the headers requested by name or expression.
        headerLines = new LinkedHashSet<VCFHeaderLine>(getSelectedHeaders(headerLines));

        // Optionally add in the reference.
        if (includeReference && getToolkit().getArguments().referenceFile != null)
            headerLines.add(new VCFHeaderLine(VCFHeader.REFERENCE_KEY, FilenameUtils.getBaseName(getToolkit().getArguments().referenceFile.getName())));

        // Optionally add in the intervals.
        if (includeIntervals && getToolkit().getArguments().intervals != null) {
            for (IntervalBinding<Feature> intervalBinding : getToolkit().getArguments().intervals) {
                String source = intervalBinding.getSource();
                if (source == null)
                    continue;
                File file = new File(source);
                if (file.exists()) {
                    headerLines.add(new VCFHeaderLine(VCFHeader.INTERVALS_KEY, FilenameUtils.getBaseName(file.getName())));
                } else {
                    headerLines.add(new VCFHeaderLine(VCFHeader.INTERVALS_KEY, source));
                }
            }
        }

        TreeSet<String> vcfSamples = new TreeSet<String>(SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE));
        VCFHeader vcfHeader = new VCFHeader(headerLines, vcfSamples);
        vcfHeader.setWriteEngineHeaders(includeEngineHeaders);
        vcfWriter.writeHeader(vcfHeader);
    }

    private Set<VCFHeaderLine> getSelectedHeaders(Set<VCFHeaderLine> headerLines) {
        Set<VCFHeaderLine> selectedHeaders = new TreeSet<VCFHeaderLine>();
        if (headerNames == null && headerExpressions == null) {
            // Include everything if nothing was explicitly included.
            selectedHeaders.addAll(headerLines);
        } else {
            // Only include the selected headers.
            if (headerNames != null)
                selectedHeaders.addAll(ListFileUtils.includeMatching(headerLines, headerKey, headerNames, true));
            if (headerExpressions != null)
                selectedHeaders.addAll(ListFileUtils.includeMatching(headerLines, headerKey, headerExpressions, false));
        }

        // Remove any excluded headers.
        if (XLheaderNames != null)
            selectedHeaders = ListFileUtils.excludeMatching(selectedHeaders, headerKey, XLheaderNames, true);
        return selectedHeaders;
    }

    /**
     * Pass through the VC record
     *
     * @param tracker the ROD tracker
     * @param ref     reference information
     * @param context alignment info
     * @return number of records processed
     */
    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        int count = 0;
        if (tracker != null) {
            Collection<VariantContext> vcs = tracker.getValues(variantCollection.variants, context.getLocation());
            if (vcs != null) {
                for (VariantContext vc : vcs) {
                    vcfWriter.add(vc);
                    count++;
                }
            }
        }
        return count;
    }

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    @Override
    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;
    }

    @Override
    public void onTraversalDone(Integer result) {
        logger.info(result + " records processed.");
    }
}
