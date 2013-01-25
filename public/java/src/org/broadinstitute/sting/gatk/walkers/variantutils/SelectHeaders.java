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

import org.apache.commons.io.FilenameUtils;
import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.interval.IntervalMergingRule;
import org.broadinstitute.sting.utils.interval.IntervalSetRule;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.variant.vcf.*;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.sting.utils.text.ListFileUtils;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextUtils;

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
@DocumentedGATKFeature( groupName = "Variant Evaluation and Manipulation Tools", extraDocs = {CommandLineGATK.class} )
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

        Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), logger);

        headerLines.add(new VCFHeaderLine(VCFHeader.SOURCE_KEY, "SelectHeaders"));

        // Select only the headers requested by name or expression.
        headerLines = new LinkedHashSet<VCFHeaderLine>(getSelectedHeaders(headerLines));

        // Optionally add in the intervals.
        if (includeIntervals) {
            IntervalArgumentCollection intervalArguments = getToolkit().getArguments().intervalArguments;
            if (intervalArguments.intervals != null) {
                for (IntervalBinding<Feature> intervalBinding : intervalArguments.intervals) {
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

            if (intervalArguments.excludeIntervals != null) {
                for (IntervalBinding<Feature> intervalBinding : intervalArguments.excludeIntervals) {
                    String source = intervalBinding.getSource();
                    if (source == null)
                        continue;
                    File file = new File(source);
                    if (file.exists()) {
                        headerLines.add(new VCFHeaderLine(VCFHeader.EXCLUDE_INTERVALS_KEY, FilenameUtils.getBaseName(file.getName())));
                    } else {
                        headerLines.add(new VCFHeaderLine(VCFHeader.EXCLUDE_INTERVALS_KEY, source));
                    }
                }
            }

            if (intervalArguments.intervalMerging != IntervalMergingRule.ALL) {
                headerLines.add(new VCFHeaderLine(VCFHeader.INTERVAL_MERGING_KEY, String.valueOf(intervalArguments.intervalMerging)));
            }

            if (intervalArguments.intervalSetRule != IntervalSetRule.UNION) {
                headerLines.add(new VCFHeaderLine(VCFHeader.INTERVAL_SET_RULE_KEY, String.valueOf(intervalArguments.intervalSetRule)));
            }

            if (intervalArguments.intervalPadding != 0) {
                headerLines.add(new VCFHeaderLine(VCFHeader.INTERVAL_PADDING_KEY, String.valueOf(intervalArguments.intervalPadding)));
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

        // always include the contig lines
        selectedHeaders = VCFUtils.withUpdatedContigsAsLines(selectedHeaders, getToolkit().getArguments().referenceFile, getToolkit().getMasterSequenceDictionary(), true);
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
