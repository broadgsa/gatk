/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils.haplotypeBAMWriter;

import net.sf.samtools.SAMFileWriter;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.haplotype.Haplotype;
import org.broadinstitute.sting.utils.genotyper.MostLikelyAllele;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.variant.variantcontext.Allele;

import java.util.*;

/**
 * Writes a BAM containing just the reads in stratifiedReadMap aligned to their
 * most likely haplotype among all of the called haplotypes.
 *
 * Primarily useful for users of the HaplotypeCaller who want to better understand the
 * support of their calls w.r.t. the reads.
 *
 * User: depristo
 * Date: 2/22/13
 * Time: 1:50 PM
 */
class CalledHaplotypeBAMWriter extends HaplotypeBAMWriter {
    public CalledHaplotypeBAMWriter(final SAMFileWriter bamWriter) {
        super(bamWriter);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void writeReadsAlignedToHaplotypes(final List<Haplotype> haplotypes,
                                              final GenomeLoc paddedReferenceLoc,
                                              final List<Haplotype> bestHaplotypes,
                                              final Set<Haplotype> calledHaplotypes,
                                              final Map<String, PerReadAlleleLikelihoodMap> stratifiedReadMap) {
        if ( calledHaplotypes.isEmpty() ) // only write out called haplotypes
            return;

        writeHaplotypesAsReads(calledHaplotypes, calledHaplotypes, paddedReferenceLoc);

        // we need to remap the Alleles back to the Haplotypes; inefficient but unfortunately this is a requirement currently
        final Map<Allele, Haplotype> alleleToHaplotypeMap = new HashMap<Allele, Haplotype>(haplotypes.size());
        for ( final Haplotype haplotype : calledHaplotypes ) {
            alleleToHaplotypeMap.put(Allele.create(haplotype.getBases()), haplotype);
        }

        // the set of all alleles that were actually called
        final Set<Allele> allelesOfCalledHaplotypes = alleleToHaplotypeMap.keySet();

        // next, output the interesting reads for each sample aligned against one of the called haplotypes
        for ( final PerReadAlleleLikelihoodMap readAlleleLikelihoodMap : stratifiedReadMap.values() ) {
            for ( Map.Entry<GATKSAMRecord, Map<Allele, Double>> entry : readAlleleLikelihoodMap.getLikelihoodReadMap().entrySet() ) {
                if ( entry.getKey().getMappingQuality() > 0 ) {
                    final MostLikelyAllele bestAllele = PerReadAlleleLikelihoodMap.getMostLikelyAllele(entry.getValue(), allelesOfCalledHaplotypes);
                    writeReadAgainstHaplotype(entry.getKey(), alleleToHaplotypeMap.get(bestAllele.getMostLikelyAllele()), paddedReferenceLoc.getStart());
                }
            }
        }
    }
}
