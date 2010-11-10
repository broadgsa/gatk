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

package org.broadinstitute.sting.gatk.walkers.variantutils;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.SortingVCFWriter;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFWriter;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.util.*;

/**
 * Left-aligns indels from a variants file.
 */
@Reference(window=@Window(start=-200,stop=200))
@Requires(value={},referenceMetaData=@RMD(name="variant", type=VariantContext.class))
public class LeftAlignVariants extends RodWalker<Integer, Integer> {

    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter baseWriter = null;

    private SortingVCFWriter writer;

    public void initialize() {
        Set<String> samples = SampleUtils.getSampleListWithVCFHeader(getToolkit(), Arrays.asList("variant"));
        Map<String, VCFHeader> vcfHeaders = VCFUtils.getVCFHeadersFromRods(getToolkit(), Arrays.asList("variant"));

        Set<VCFHeaderLine> headerLines = vcfHeaders.get("variant").getMetaData();
        baseWriter.writeHeader(new VCFHeader(headerLines, samples));

        writer = new SortingVCFWriter(baseWriter, 200);
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        Collection<VariantContext> VCs = tracker.getVariantContexts(ref, "variant", null, context.getLocation(), true, false);

        int changedSites = 0;
        for ( VariantContext vc : VCs )
            changedSites += alignAndWrite(vc, ref);

        return changedSites;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {
        System.out.println(result + " variants were aligned");
    }


    private int alignAndWrite(VariantContext vc, final ReferenceContext ref) {
        if ( vc.isBiallelic() && vc.isIndel() )
            return writeLeftAlignedIndel(vc, ref);
        else {
            writer.add(vc, ref.getBase());
            return 0;
        }
    }

    private int writeLeftAlignedIndel(final VariantContext vc, final ReferenceContext ref) {
        final byte[] refSeq = ref.getBases();

        // get the indel length
        int indelLength;
        if ( vc.isDeletion() )
            indelLength = vc.getReference().length();
        else
            indelLength = vc.getAlternateAllele(0).length();

        // create an indel haplotype
        int originalIndex = (int)(ref.getLocus().getStart() - ref.getWindow().getStart()) + 1;
        final byte[] originalIndel = makeHaplotype(vc, refSeq, originalIndex, indelLength);

        // create a CIGAR string to represent the event
        ArrayList<CigarElement> elements = new ArrayList<CigarElement>();
        elements.add(new CigarElement(originalIndex, CigarOperator.M));
        elements.add(new CigarElement(indelLength, vc.isDeletion() ? CigarOperator.D : CigarOperator.I));
        elements.add(new CigarElement(refSeq.length - originalIndex, CigarOperator.M));
        Cigar originalCigar = new Cigar(elements);

        // left align the CIGAR
        Cigar newCigar = AlignmentUtils.leftAlignIndel(originalCigar, refSeq, originalIndel, 0, 0);

        // update if necessary and write
        if ( !newCigar.equals(originalCigar) && newCigar.numCigarElements() > 1 ) {
            int difference = originalIndex - newCigar.getCigarElement(0).getLength();
            GenomeLoc newLoc = getToolkit().getGenomeLocParser().createPotentiallyInvalidGenomeLoc(vc.getChr(), vc.getStart()-difference, vc.getEnd()-difference);
            //System.out.println("Moving record from " + vc.getChr()+":"+vc.getStart() + " to " + newLoc);
            VariantContext newVC = VariantContextUtils.modifyLocation(vc, newLoc);

            int indelIndex = originalIndex-difference;
            byte[] newBases = new byte[indelLength];
            System.arraycopy(refSeq, indelIndex, newBases, 0, indelLength);
            Allele newAllele = Allele.create(newBases, vc.isDeletion());
            newVC = updateAllele(newVC, newAllele);

            writer.add(newVC, refSeq[indelIndex-1]);
            return 1;
        } else {
            writer.add(vc, ref.getBase());
            return 0;
        }
    }

    private static byte[] makeHaplotype(VariantContext vc, byte[] ref, int indexOfRef, int indelLength) {
        byte[] hap = new byte[ref.length + (indelLength * (vc.isDeletion() ? -1 : 1))];

        // add the bases before the indel
        System.arraycopy(ref, 0, hap, 0, indexOfRef);
        int currentPos = indexOfRef;

        // take care of the indel
        if ( vc.isDeletion() ) {
            indexOfRef += indelLength;
        } else {
            System.arraycopy(vc.getAlternateAllele(0).getBases(), 0, hap, currentPos, indelLength);
            currentPos += indelLength;
        }

        // add the bases after the indel
        System.arraycopy(ref, indexOfRef, hap, currentPos, ref.length - indexOfRef);

        return hap;
    }

    public static VariantContext updateAllele(VariantContext vc, Allele newAllele) {
        // create a mapping from original allele to new allele
        HashMap<Allele, Allele> alleleMap = new HashMap<Allele, Allele>(vc.getAlleles().size());
        if ( newAllele.isReference() ) {
            alleleMap.put(vc.getReference(), newAllele);
            alleleMap.put(vc.getAlternateAllele(0), vc.getAlternateAllele(0));
        } else {
            alleleMap.put(vc.getReference(), vc.getReference());
            alleleMap.put(vc.getAlternateAllele(0), newAllele);
        }

        // create new Genotype objects
        Map<String, Genotype> newGenotypes = new HashMap<String, Genotype>(vc.getNSamples());
        for ( Map.Entry<String, Genotype> genotype : vc.getGenotypes().entrySet() ) {
            List<Allele> newAlleles = new ArrayList<Allele>();
            for ( Allele allele : genotype.getValue().getAlleles() ) {
                Allele newA = alleleMap.get(allele);
                if ( newA == null )
                    newA = Allele.NO_CALL;
                newAlleles.add(newA);
            }
            newGenotypes.put(genotype.getKey(), Genotype.modifyAlleles(genotype.getValue(), newAlleles));
        }

        return new VariantContext(vc.getSource(), vc.getChr(), vc.getStart(), vc.getEnd(), alleleMap.values(), newGenotypes, vc.getNegLog10PError(), vc.filtersWereApplied() ? vc.getFilters() : null, vc.getAttributes());
    }
}
