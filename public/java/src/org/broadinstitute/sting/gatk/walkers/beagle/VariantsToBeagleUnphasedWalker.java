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

package org.broadinstitute.sting.gatk.walkers.beagle;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Set;

/**
 * Produces an input file to Beagle imputation engine, listing unphased, hard-called genotypes for a single sample
 * in input variant file.  Will additional hold back a fraction of the sites for evaluation, marking the
 * genotypes at that sites as missing, and writing the truth of these sites to a second VCF file
 */
public class VariantsToBeagleUnphasedWalker extends RodWalker<Integer, Integer> {
    @Input(fullName="variants", shortName = "V", doc="Input VCF file", required=true)
    public RodBinding<VariantContext> variants;

    @Output(doc="File to which BEAGLE unphased genotypes should be written",required=true)
    protected PrintStream  beagleWriter = null;

    @Argument(fullName = "bootstrap_fraction", shortName = "bs", doc = "Proportion of records to be used in bootstrap set", required = false)
    public double bootstrap = 0.0;

    @Argument(fullName = "bootstrap_vcf",shortName = "bsvcf", doc = "Output a VCF with the records used for bootstrapping filtered out", required = false)
    VCFWriter bootstrapVCFOutput = null;

    @Argument(fullName = "missing", shortName = "missing", doc = "String to identify missing data in beagle output", required = false)
    public String MISSING = "?";

    private Set<String> samples = null;
    private int bootstrapSetSize = 0;
    private int testSetSize = 0;

    public void initialize() {
        samples = SampleUtils.getSampleListWithVCFHeader(getToolkit(), Arrays.asList(variants.getName()));

        beagleWriter.print("I marker alleleA alleleB");
        for ( String sample : samples )
            beagleWriter.print(String.format(" %s %s", sample, sample));

        beagleWriter.println();

        if ( bootstrap < 0.0 | bootstrap > 1.0 )
            throw new UserException.BadArgumentValue("bootstrap", "Bootstrap value must be fraction between 0 and 1");

        if ( bootstrapVCFOutput != null ) {
            Set<VCFHeaderLine> hInfo = VCFUtils.getHeaderFields(getToolkit());
            bootstrapVCFOutput.writeHeader(new VCFHeader(hInfo, SampleUtils.getUniqueSamplesFromRods(getToolkit())));
        }
    }

    /**
     * Iterate over each site, emitting the BEAGLE unphased genotypes file format
     * @param tracker
     * @param ref
     * @param context
     * @return
     */
    public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        if( tracker != null ) {
            GenomeLoc loc = context.getLocation();
            VariantContext vc = tracker.getFirstValue(variants, loc);

            if ( ProduceBeagleInputWalker.canBeOutputToBeagle(vc) ) {
                // do we want to hold back this site?
                boolean makeMissing = dropSite(vc);

                // if we are holding it back and we are writing a bootstrap VCF, write it out
                if ( makeMissing && bootstrapVCFOutput != null ) {
                    bootstrapVCFOutput.add(vc);
                }

                // regardless, all sites are written to the unphased genotypes file, marked as missing if appropriate
                writeUnphasedBeagleOutput(vc, makeMissing);
            }
        }

        return 0;
    }

    /**
     * Do we want to hold back this site for bootstrap?  Considers the bootstrap fraction member variable
     *
     * @param vc
     * @return
     */
    public boolean dropSite(VariantContext vc) {
        if ( (bootstrapSetSize+1.0)/(1.0+bootstrapSetSize+testSetSize) <= bootstrap ) {
            bootstrapSetSize++;
            return true;
        } else {
            testSetSize++;
            return false;
        }
    }

    public void writeUnphasedBeagleOutput(VariantContext vc, boolean makeMissing) {
        GenomeLoc currentLoc = VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(),vc);
        StringBuffer beagleOut = new StringBuffer();

        String marker = String.format("%s:%d ",currentLoc.getContig(), currentLoc.getStart());
        beagleOut.append("M ").append(marker);

        // write out the alleles at this site
        for ( Allele allele : vc.getAlleles() ) {
            beagleOut.append(allele.isNoCall() || allele.isNull() ? "-" : allele.getBaseString()).append(" ");
        }

        // write out sample level genotypes
        for ( String sample : samples ) {
            Genotype genotype = vc.getGenotype(sample);
            if ( ! makeMissing && genotype.isCalled() ) {
                addAlleles(beagleOut, genotype);
            } else {
                addAlleles(beagleOut, MISSING, MISSING);
            }
        }

        beagleWriter.println(beagleOut.toString());
    }

    private void addAlleles(StringBuffer buf, Genotype g) {
        addAlleles(buf, g.getAllele(0).getBaseString(), g.getAllele(1).getBaseString());

    }

    private void addAlleles(StringBuffer buf, String a, String b) {
        buf.append(a).append(" ").append(b);
    }

    public Integer reduceInit() { return 0; }
    public Integer reduce( Integer value, Integer sum ) { return value + sum; }

    public void onTraversalDone( Integer includedSites ) {
        logger.info("Sites included in beagle genotypes file : " + includedSites);
    }
}
