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

package org.broadinstitute.gatk.tools.walkers.varianteval.stratifications;

import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.gatk.tools.walkers.varianteval.evaluators.VariantSummary;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;

import java.util.*;

/**
 * Stratifies the eval RODs by the allele count of the alternate allele
 *
 * Looks first at the MLEAC value in the INFO field, and uses that value if present.
 * If not present, it then looks for the AC value in the INFO field.  If both are absent,
 * it computes the AC from the genotypes themselves.  If no AC can be computed, 0 is used.
 */
public class AlleleCount extends VariantStratifier {
    int nchrom;

    @Override
    public void initialize() {
        // we can only work with a single eval VCF, and it must have genotypes
        if ( getVariantEvalWalker().getEvals().size() != 1 && !getVariantEvalWalker().mergeEvals )
            throw new UserException.BadArgumentValue("AlleleCount", "AlleleCount stratification only works with a single eval vcf");

        // There are ploidy x n sample chromosomes
        // TODO -- generalize to handle multiple ploidy
        nchrom = getVariantEvalWalker().getNumberOfSamplesForEvaluation() * getVariantEvalWalker().getSamplePloidy();
        if ( nchrom < 2 )
            throw new UserException.BadArgumentValue("AlleleCount", "AlleleCount stratification requires an eval vcf with at least one sample");

        // create an array containing each of the allele counts
        for( int ac = 0; ac <= nchrom; ac++ ) {
            states.add(ac);
        }

        getVariantEvalWalker().getLogger().info("AlleleCount using " + nchrom + " chromosomes");
    }

    public List<Object> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String familyName) {
        if (eval != null) {
            int AC = 0; // by default, the site is considered monomorphic

            try {
                if ( eval.isBiallelic() ) {
                    if ( eval.hasAttribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY) ) {
                        // the MLEAC is allowed to be larger than the AN (e.g. in the case of all PLs being 0, the GT is ./. but the exact model may arbitrarily choose an AC>1)
                        AC = Math.min(eval.getAttributeAsInt(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, 0), nchrom);
                    } else if ( eval.hasAttribute(VCFConstants.ALLELE_COUNT_KEY) ) {
                        AC = eval.getAttributeAsInt(VCFConstants.ALLELE_COUNT_KEY, 0);
                    }
                }
            } catch ( ClassCastException e ) {
                // protect ourselves from bad inputs
                // TODO -- fully decode VC
            }

            if ( AC == 0 && eval.isVariant() ) {
                // fall back to the direct calculation
                for (Allele allele : eval.getAlternateAlleles())
                    AC = Math.max(AC, eval.getCalledChrCount(allele));
            }

            // make sure that the AC isn't invalid
            if ( AC > nchrom )
                throw new UserException.MalformedVCF(String.format("The AC value (%d) at position %s:%d " +
                        "is larger than the number of chromosomes over all samples (%d)", AC,
                        eval.getChr(), eval.getStart(), nchrom));

            return Collections.singletonList((Object) AC);
        } else {
            return Collections.emptyList();
        }
    }

    @Override
    public Set<Class<? extends VariantEvaluator>> getIncompatibleEvaluators() {
        return new HashSet<Class<? extends VariantEvaluator>>(Arrays.asList(VariantSummary.class));
    }

    @Override
    public String getFormat() {
        return "%d";
    }
}
