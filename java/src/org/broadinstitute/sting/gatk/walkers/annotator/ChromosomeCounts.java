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

package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broad.tribble.vcf.VCFRecord;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.*;

import java.util.*;


public class ChromosomeCounts implements InfoFieldAnnotation, StandardAnnotation {

    private String[] keyNames = { VCFRecord.ALLELE_NUMBER_KEY, VCFRecord.ALLELE_COUNT_KEY, VCFRecord.ALLELE_FREQUENCY_KEY };
    private VCFInfoHeaderLine[] descriptions = { new VCFInfoHeaderLine(VCFRecord.ALLELE_FREQUENCY_KEY, -1, VCFInfoHeaderLine.INFO_TYPE.Float, "Allele Frequency"),
            new VCFInfoHeaderLine(VCFRecord.ALLELE_COUNT_KEY, -1, VCFInfoHeaderLine.INFO_TYPE.Integer, "Allele count in genotypes, for each ALT allele, in the same order as listed"),
            new VCFInfoHeaderLine(VCFRecord.ALLELE_NUMBER_KEY, 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "Total number of alleles in called genotypes") };

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {

        if ( vc.getChromosomeCount() == 0 )
            return null;
        
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(VCFRecord.ALLELE_NUMBER_KEY, vc.getChromosomeCount());

        if ( vc.getAlternateAlleles().size() > 0 ) {
            ArrayList<Double> alleleFreqs = new ArrayList<Double>();
            ArrayList<Integer> alleleCounts = new ArrayList<Integer>();
            for ( Allele allele : vc.getAlternateAlleles() ) {
                alleleCounts.add(vc.getChromosomeCount(allele));
                alleleFreqs.add((double)vc.getChromosomeCount(allele) / (double)vc.getChromosomeCount());
            }

            map.put(VCFRecord.ALLELE_COUNT_KEY, alleleCounts);
            map.put(VCFRecord.ALLELE_FREQUENCY_KEY, alleleFreqs);
        }
        
        return map;
    }

    public List<String> getKeyNames() {
        return Arrays.asList(keyNames);
    }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(descriptions); }
}