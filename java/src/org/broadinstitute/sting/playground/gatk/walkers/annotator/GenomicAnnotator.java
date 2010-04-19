/*
 * Copyright (c) 2010 The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the ”Software”), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED ”AS IS”, WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.gatk.walkers.annotator;

import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.walkers.Allows;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.genotype.vcf.VCFHeader;
import org.broadinstitute.sting.utils.genotype.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.genotype.vcf.VCFUtils;
import org.broadinstitute.sting.utils.genotype.vcf.VCFWriter;


/**
 * Annotates variant calls with information from user-specified tabular files.
 *
 * For details, see:  http://www.broadinstitute.org/gsa/wiki/index.php/GenomicAnnotator
 */
//@Requires(value={DataSource.READS, DataSource.REFERENCE},referenceMetaData=@RMD(name="variant",type=VariationRod.class))
@Allows(value={DataSource.READS, DataSource.REFERENCE})
@Reference(window=@Window(start=-50,stop=50))
@By(DataSource.REFERENCE)
public class GenomicAnnotator extends RodWalker<Integer, Integer> {
    @Argument(fullName="vcfOutput", shortName="vcf", doc="VCF file to which all variants should be written with annotations", required=true)
    protected File VCF_OUT;

    @Argument(fullName="sampleName", shortName="sample", doc="The sample (NA-ID) corresponding to the variant input (for non-VCF input only)", required=false)
    protected String sampleName = null;

    @Argument(fullName="select", shortName="s", doc="Optionally specifies which subset of columns from which -B inputs should be used for annotations. For example, -B mydbsnp,AnnotatorInputTable,/path/to/mydbsnp.txt -B mytable,AnnotatorInputTable,/path/mytable.txt -s mydbsnp.avHet,mydbsnp.name,mytable.column3 will cause annotations to only be generated from the 3 columns specified using -s.", required=false)
    protected String[] COLUMNS = {};

    @Argument(fullName="oneToMany", shortName="m", doc="If more than one record from the same file matches a particular locus (for example, multiple dbSNP records with the same position), create multiple entries in the ouptut VCF file - one for each match. If a particular tabular file has J matches, and another tabular file has K matches for a given locus, then J*K output VCF records will be generated - one for each pair of K, J.   If this flag is not provided, the multiple records are still generated, but they are stored in the INFO field of a single output VCF record, with their annotation keys differentiated by appending '_i' with i varying from 1 to K*J. ", required=false)
    protected Boolean ONE_TO_MANY = false;

    private VCFWriter vcfWriter;

    private HashMap<String, String> nonVCFsampleName = new HashMap<String, String>();

    private VariantAnnotatorEngine engine;


    /**
     * Prepare the output file and the list of available features.
     */
    public void initialize() {

        // get the list of all sample names from the various VCF input rods
        TreeSet<String> samples = new TreeSet<String>();
        SampleUtils.getUniquifiedSamplesFromRods(getToolkit(), samples, new HashMap<Pair<String, String>, String>());

        // add the non-VCF sample from the command-line, if applicable
        if ( sampleName != null  ) {
            nonVCFsampleName.put(sampleName.toUpperCase(), "variant");
            samples.add(sampleName.toUpperCase());
        }

        // if there are no valid samples, warn the user
        if ( samples.size() == 0 ) {
            logger.warn("There are no samples input at all; use the --sampleName argument to specify one if desired.");
        }

        engine = new VariantAnnotatorEngine(getToolkit(), new String[] { }, new String[] { "GenomicAnnotation" });

        engine.setOneToMany( Boolean.TRUE.equals( ONE_TO_MANY ) );
        engine.setRequestedColumns(COLUMNS);

        // setup the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("source", "Annotator"));
        hInfo.add(new VCFHeaderLine("annotatorReference", getToolkit().getArguments().referenceFile.getName()));
        hInfo.addAll(engine.getVCFAnnotationDescriptions());

        vcfWriter = new VCFWriter(VCF_OUT);
        VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);
    }

    /**
     * Initialize the number of loci processed to zero.
     *
     * @return 0
     */
    public Integer reduceInit() { return 0; }


    /**
     * We want reads that span deletions
     *
     * @return true
     */
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    /**
     * For each site of interest, annotate based on the requested annotation types
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return 1 if the locus was successfully processed, 0 if otherwise
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        List<Object> rods = tracker.getReferenceMetaData("variant");
        // ignore places where we don't have a variant
        if ( rods.size() == 0 )
            return 0;

        Object variant = rods.get(0);
        if( ref.getBase() == 'N') {
            return 0; //TODO Currently, VariantContextAdaptors.toVCF(annotatedVC, ref.getBase()) fails when base is 'N'. is this right?
        }

        VariantContext vc = VariantContextAdaptors.toVariantContext("variant", variant);
        if ( vc == null )
            return 0;

        // if the reference base is not ambiguous, we can annotate
        Collection<VariantContext> annotatedVCs = Arrays.asList(vc);
        if ( BaseUtils.simpleBaseToBaseIndex(ref.getBase()) != -1 ) {
            Map<String, StratifiedAlignmentContext> stratifiedContexts = StratifiedAlignmentContext.splitContextBySample(context.getBasePileup());
            if ( stratifiedContexts != null ) {
                annotatedVCs = engine.annotateContext(tracker, ref, stratifiedContexts, vc);
            }
        }

        if ( variant instanceof RodVCF ) { //TODO is this if-statement necessary?
            for(VariantContext annotatedVC : annotatedVCs ) {
                    vcfWriter.addRecord(VariantContextAdaptors.toVCF(annotatedVC, ref.getBase()));
            }
        }

        return 1;
    }

    /**
     * Increment the number of loci processed.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return the new number of loci processed.
     */
    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    /**
     * Tell the user the number of loci processed and close out the new variants file.
     *
     * @param result  the number of loci seen.
     */
    public void onTraversalDone(Integer result) {
        out.printf("Processed %d loci.\n", result);

        vcfWriter.close();
    }
}

