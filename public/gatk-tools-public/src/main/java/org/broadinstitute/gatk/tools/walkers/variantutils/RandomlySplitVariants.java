/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.variantutils;

import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterFactory;

import java.io.File;
import java.util.*;

/**
 * Randomly split variants into different sets
 *
 * <p>This tool takes a VCF file, randomly splits variants into different sets, and writes the
 * results to separate files. By default the tool splits the input into two new sets, but it can be made to output
 * more than two separate call sets.</p>
 *
 * <h3>Input</h3>
 * <p>
 * A variant call set to split.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * The new callsets.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T RandomlySplitVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -o1 output_1.vcf \
 *   -o2 output_2.vcf
 * </pre>
 *
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
public class RandomlySplitVariants extends RodWalker<Integer, Integer> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Output(fullName="out1", shortName="o1", doc="File #1 to which variants should be written", required=false, exclusiveOf = "splitToManyFiles")
    protected VariantContextWriter vcfWriter1 = null;

    @Output(fullName="out2", shortName="o2", doc="File #2 to which variants should be written", required=false, exclusiveOf = "splitToManyFiles")
    // there's a reported bug in the GATK where we can't have 2 @Output writers
    protected File file2 = null;
    protected VariantContextWriter vcfWriter2 = null;

    @Argument(fullName="fractionToOut1", shortName="fraction", doc="Fraction of records to be placed in out1 (must be 0 >= fraction <= 1); all other records are placed in out2", required=false)
    protected double fraction = 0.5;

    @Argument(fullName="splitToManyFiles", shortName = "splitToMany", doc="split (with uniform distribution) to more than 2 files. numOfFiles and baseOutputName parameters are required", required = false)
    protected boolean splitToMany = false;

    @Argument(fullName = "numOfOutputVCFFiles", shortName = "N", doc = "number of output VCF files. Only works with SplitToMany = true", required = false, maxRecommendedValue = 20, minValue = 2)
    protected int numOfFiles = -1;

    @Argument(fullName = "prefixForAllOutputFileNames", shortName = "baseOutputName", doc = "the name of the output VCF file will be: <baseOutputName>.split.<number>.vcf. Required with SplitToMany option", required = false)
    protected String baseFileName = null;

    private VariantContextWriter[] writers = null;

    /**
     * Set up the VCF writer, the sample expressions and regexs, and the JEXL matcher
     */
    public void initialize() {
        if ( fraction < 0.0 || fraction > 1.0 )
            throw new UserException.BadArgumentValue("fractionToOut1", "this value needs to be a number between 0 and 1");

        if (splitToMany){
            if (numOfFiles < 2)
                throw new UserException.BadArgumentValue("numOfFiles", "this value must be greater than 2 when using the splitToMany option");
            if (baseFileName == null)
                throw new UserException.BadArgumentValue("baseFileName", "this value cannot be null (unprovided) when using the splitToMany option");
        }
        else{
            if(vcfWriter1 == null || vcfWriter2 == null)
                throw new UserException.BadArgumentValue("out1 or out2", "this value cannot be null (unprovided) unless you are using the splitToMany option");
        }

        // setup the header info
        final List<String> inputNames = Arrays.asList(variantCollection.variants.getName());
        final Set<String> samples = SampleUtils.getUniqueSamplesFromRods(getToolkit(), inputNames);
        final Set<VCFHeaderLine> hInfo = new HashSet<>();
        hInfo.addAll(GATKVCFUtils.getHeaderFields(getToolkit(), inputNames));


        if(splitToMany){
            writers = new VariantContextWriter[numOfFiles];
            for(int i = 0; i<writers.length; i++){
                writers[i] = VariantContextWriterFactory.create(new File(baseFileName+".split."+i+".vcf"), getMasterSequenceDictionary());
                writers[i].writeHeader(new VCFHeader(hInfo,samples));
            }

        }
        else {
            vcfWriter1.writeHeader(new VCFHeader(hInfo, samples));
            vcfWriter2 = VariantContextWriterFactory.create(file2, getMasterSequenceDictionary());
            vcfWriter2.writeHeader(new VCFHeader(hInfo, samples));
        }
    }

    /**
     * Subset VC record if necessary and emit the modified record (provided it satisfies criteria for printing)
     *
     * @param  tracker   the ROD tracker
     * @param  ref       reference information
     * @param  context   alignment info
     * @return 1 if the record was printed to the output file, 0 if otherwise
     */
    public Integer map(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {
        if ( tracker == null )
            return 0;

        final Collection<VariantContext> vcs = tracker.getValues(variantCollection.variants, context.getLocation());
        for ( final VariantContext vc : vcs ) {
            final double random = Utils.getRandomGenerator().nextDouble();
            if(splitToMany){
                final int index = (int)(numOfFiles * random);
                writers[index].add(vc);
            }
            else{
                if ( random < fraction )
                    vcfWriter1.add(vc);
                else
                    vcfWriter2.add(vc);
            }
        }

        return 1;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(final Integer value, final Integer sum) { return value + sum; }

    public void onTraversalDone(final Integer result) {
        logger.info(result + " records processed.");
        if(splitToMany)
            for(final VariantContextWriter writer: writers)
                writer.close();
        else
            vcfWriter2.close();
    }
}
