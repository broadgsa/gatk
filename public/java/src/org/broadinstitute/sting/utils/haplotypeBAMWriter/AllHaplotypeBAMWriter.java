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

import net.sf.samtools.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.haplotype.Haplotype;
import org.broadinstitute.sting.utils.genotyper.MostLikelyAllele;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.variant.variantcontext.Allele;

import java.util.*;

/**
 * A haplotype bam writer that writes out all haplotypes as reads and then
 * the alignment of reach read to its best match among the best haplotypes.
 *
 * Primarily useful for people working on the HaplotypeCaller method itself
 *
 * User: depristo
 * Date: 2/22/13
 * Time: 1:50 PM
 */
class AllHaplotypeBAMWriter extends HaplotypeBAMWriter {
    public AllHaplotypeBAMWriter(final SAMFileWriter bamWriter) {
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
        writeHaplotypesAsReads(haplotypes, new HashSet<Haplotype>(bestHaplotypes), paddedReferenceLoc);

        // we need to remap the Alleles back to the Haplotypes; inefficient but unfortunately this is a requirement currently
        final Map<Allele, Haplotype> alleleToHaplotypeMap = new HashMap<Allele, Haplotype>(haplotypes.size());
        for ( final Haplotype haplotype : haplotypes )
            alleleToHaplotypeMap.put(Allele.create(haplotype.getBases()), haplotype);

        // next, output the interesting reads for each sample aligned against the appropriate haplotype
        for ( final PerReadAlleleLikelihoodMap readAlleleLikelihoodMap : stratifiedReadMap.values() ) {
            for ( Map.Entry<GATKSAMRecord, Map<Allele, Double>> entry : readAlleleLikelihoodMap.getLikelihoodReadMap().entrySet() ) {
                final MostLikelyAllele bestAllele = PerReadAlleleLikelihoodMap.getMostLikelyAllele(entry.getValue());
                writeReadAgainstHaplotype(entry.getKey(), alleleToHaplotypeMap.get(bestAllele.getMostLikelyAllele()), paddedReferenceLoc.getStart());
            }
        }
    }
}
