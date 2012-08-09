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

import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class UniformSamplingFrequencySelector extends FrequencyModeSelector {
    private ArrayList<GenomeEvent> binnedEventArray;

    public UniformSamplingFrequencySelector(GenomeLocParser parser) {
        super(parser);
        binnedEventArray = new ArrayList<GenomeEvent>();

    }

    public void logCurrentSiteData(VariantContext vc, boolean selectedInTargetSamples, boolean IGNORE_GENOTYPES, boolean IGNORE_POLYMORPHIC) {
        HashMap<String, Object> attributes = new HashMap<String, Object>();


        if (vc.hasGenotypes() && !IGNORE_GENOTYPES) {
            // recompute AF,AC,AN based on genotypes:
            VariantContextUtils.calculateChromosomeCounts(vc, attributes, false);
            if (! selectedInTargetSamples && !IGNORE_POLYMORPHIC)
                return;
        } else  {
            if (!IGNORE_POLYMORPHIC) {
                if (vc.getAttributes().containsKey(VCFConstants.ALLELE_COUNT_KEY))
                {
                    int ac = vc.getAttributeAsInt(VCFConstants.ALLELE_COUNT_KEY, 0);
                    if (ac == 0) return; // site not polymorphic
                }
                else
                    // no allele count field in VC
                    return;
            }
        }
        // create bare-bones event and log in corresponding bin
        // attributes contains AC,AF,AN pulled from original vc, and we keep them here and log in output file for bookkeeping purposes
        GenomeEvent event = new GenomeEvent(parser, vc.getChr(), vc.getStart(), vc.getEnd(),vc.getAlleles(), attributes);
        binnedEventArray.add(event);

    }

    public ArrayList<VariantContext> selectValidationSites(int numValidationSites) {

        // take randomly sitesToChoosePerBin[k] elements from each bin
        ArrayList<GenomeEvent> selectedEvents = new ArrayList<GenomeEvent>();

        selectedEvents.addAll(MathUtils.randomSubset(binnedEventArray, numValidationSites));

        Collections.sort(selectedEvents);

        // now convert to VC
        ArrayList<VariantContext> selectedSites = new ArrayList<VariantContext>();
        for (GenomeEvent event : selectedEvents)
            selectedSites.add(event.createVariantContextFromEvent());

        return selectedSites;
    }
}
