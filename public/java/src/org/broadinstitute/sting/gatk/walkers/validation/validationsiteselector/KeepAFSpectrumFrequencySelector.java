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
package org.broadinstitute.sting.gatk.walkers.validation.validationsiteselector;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;


public class KeepAFSpectrumFrequencySelector extends FrequencyModeSelector {

    private static final boolean DEBUG = true;

    private int NUM_BINS = 20;

    private int[] preSampleSelectionHistogram;
    private int numTotalSites = 0;
    private int[] postSampleSelectionHistogram;
    private int numSampleSelectedSites = 0;
    private ArrayList<GenomeEvent>[] binnedEventArray;

    public KeepAFSpectrumFrequencySelector(int numBins, GenomeLocParser parser) {
        super(parser);
        NUM_BINS = numBins;
        // initialize arrays dependent on NUM_BINS
        binnedEventArray = new ArrayList[NUM_BINS];

        for (int k=0; k < NUM_BINS; k++)
            binnedEventArray[k] = new ArrayList<GenomeEvent>();

        preSampleSelectionHistogram = new int[NUM_BINS];
        postSampleSelectionHistogram = new int[NUM_BINS];
    }

    public void logCurrentSiteData(VariantContext vc, boolean selectedInTargetSamples, boolean IGNORE_GENOTYPES, boolean IGNORE_POLYMORPHIC) {

        // this method is called for every  variant of a selected type, regardless of  whether it will be selectable or not
        // get AC,AF,AN attributes from vc
        HashMap<String, Object> attributes = new HashMap<String, Object>();
        double[] afArray = null;

        if (vc.hasGenotypes() && !IGNORE_GENOTYPES) {
            // recompute AF,AC,AN based on genotypes:
            // todo - - maybe too inefficient??
            VariantContextUtils.calculateChromosomeCounts(vc, attributes, false);
        }

        // sites-only vc or we explicitly tell to ignore genotypes; we trust the AF field if present
        if ( vc.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY) )  {
            String afo = vc.getAttributeAsString(VCFConstants.ALLELE_FREQUENCY_KEY, null);

            if (afo.contains(",")) {
                String[] afs = afo.split(",");
                afs[0] = afs[0].substring(1,afs[0].length());
                afs[afs.length-1] = afs[afs.length-1].substring(0,afs[afs.length-1].length()-1);

                afArray = new double[afs.length];

                for (int k=0; k < afArray.length; k++)
                    afArray[k] = Double.valueOf(afs[k]);
            }
            else
                afArray = new double[] {Double.valueOf(afo)};
        }


        if (afArray == null )
            return;

        double af0 = MathUtils.arrayMax(afArray);

        int binIndex = (NUM_BINS-1) - (int) Math.floor(((1.0-af0)*NUM_BINS));
        // deal with round-off issue: low-AC sites with large samples can have AF rounded down to 0.000
        if (binIndex < 0)
            binIndex = 0;
//        System.out.format("Pre:%4.4f %d\n",af0, binIndex);
        preSampleSelectionHistogram[binIndex]++;
        numTotalSites++;

        // now process VC subsetted to samples of interest
        if (! selectedInTargetSamples && !IGNORE_POLYMORPHIC)
            return;

        //System.out.format("Post:%4.4f %d\n",af0, binIndex);
        postSampleSelectionHistogram[binIndex]++;
        numSampleSelectedSites++;

        // create bare-bones event and log in corresponding bin
        // attributes contains AC,AF,AN pulled from original vc, and we keep them here and log in output file for bookkeeping purposes
        GenomeEvent event = new GenomeEvent(parser, vc.getChr(), vc.getStart(), vc.getEnd(),vc.getAlleles(), attributes);

        binnedEventArray[binIndex].add(event);

    }

    public ArrayList<VariantContext> selectValidationSites(int numValidationSites) {
        // number of sites to choose at random for each frequency bin = #desired validation sites/# total sites * #sites in original bin
        int[] sitesToChoosePerBin = new int[NUM_BINS];
        int totalSites = 0;
        for (int k=0; k < NUM_BINS; k++) {
            int sites = (int)Math.round((double)numValidationSites * preSampleSelectionHistogram[k]/ (double)numTotalSites);
            sitesToChoosePerBin[k] = sites;
            totalSites += sites;
        }

        // deal with rounding artifacts
        while (totalSites > numValidationSites) {
            // take off one from randomly selected bin
            int k= GenomeAnalysisEngine.getRandomGenerator().nextInt(NUM_BINS);
            sitesToChoosePerBin[k]--;
            totalSites--;
        }
        while (totalSites < numValidationSites) {
            // take off one from randomly selected bin
            int k= GenomeAnalysisEngine.getRandomGenerator().nextInt( NUM_BINS);
            sitesToChoosePerBin[k]++;
            totalSites++;
        }

        if (DEBUG) {
            System.out.println("sitesToChoosePerBin:");
            for (int k=0; k < NUM_BINS; k++)
                System.out.format("%d ", sitesToChoosePerBin[k]);
            System.out.println();

            System.out.println("preSampleSelectionHistogram:");
            for (int k=0; k < NUM_BINS; k++)
                System.out.format("%d ", preSampleSelectionHistogram[k]);
            System.out.println();

            System.out.println("postSampleSelectionHistogram:");
            for (int k=0; k < NUM_BINS; k++)
                System.out.format("%d ", postSampleSelectionHistogram[k]);
            System.out.println();

        }

        // take randomly sitesToChoosePerBin[k] elements from each bin
        ArrayList<GenomeEvent> selectedEvents = new ArrayList<GenomeEvent>();

        for (int k=0; k < NUM_BINS; k++) {
            selectedEvents.addAll(MathUtils.randomSubset(binnedEventArray[k], sitesToChoosePerBin[k]));
        }

        Collections.sort(selectedEvents);

        // now convert to VC
        ArrayList<VariantContext> selectedSites = new ArrayList<VariantContext>();
        for (GenomeEvent event : selectedEvents)
            selectedSites.add(event.createVariantContextFromEvent());

        return selectedSites;

    }

}
