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
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

/**
 * Left-aligns indels from a variants file.
 *
 * <p>
 * LeftAlignVariants is a tool that takes a VCF file and left-aligns any indels inside it.  The same indel can often be
 * placed at multiple positions and still represent the same haplotype.  While the standard convention with VCF is to
 * place an indel at the left-most position this doesn't always happen, so this tool can be used to left-align them.
 *
 * <h2>Input</h2>
 * <p>
 * A variant set to left-align.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A left-aligned VCF.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T LeftAlignVariants \
 *   --variant input.vcf \
 *   -o output.vcf
 * </pre>
 *
 */
@Reference(window=@Window(start=-200,stop=200))
public class LeftAlignVariants extends RodWalker<Integer, Integer> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter baseWriter = null;

    private SortingVCFWriter writer;

    public void initialize() {
        String trackName = variantCollection.variants.getName();
        Set<String> samples = SampleUtils.getSampleListWithVCFHeader(getToolkit(), Arrays.asList(trackName));
        Map<String, VCFHeader> vcfHeaders = VCFUtils.getVCFHeadersFromRods(getToolkit(), Arrays.asList(trackName));

        Set<VCFHeaderLine> headerLines = vcfHeaders.get(trackName).getMetaData();
        baseWriter.writeHeader(new VCFHeader(headerLines, samples));

        writer = new SortingVCFWriter(baseWriter, 200);
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        Collection<VariantContext> VCs = tracker.getValues(variantCollection.variants, context.getLocation());

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
        writer.close();
        System.out.println(result + " variants were aligned");
    }


    private int alignAndWrite(VariantContext vc, final ReferenceContext ref) {
        if ( vc.isBiallelic() && vc.isIndel() && !vc.isComplexIndel() )
            return writeLeftAlignedIndel(vc, ref);
        else {
            writer.add(vc);
            return 0;
        }
    }

    private int writeLeftAlignedIndel(final VariantContext vc, final ReferenceContext ref) {
        final byte[] refSeq = ref.getBases();

        // get the indel length
        int indelLength;
        if ( vc.isSimpleDeletion() )
            indelLength = vc.getReference().length();
        else
            indelLength = vc.getAlternateAllele(0).length();

        if ( indelLength > 200 ) {
            writer.add(vc);
            return 0;
        }

        // create an indel haplotype
        int originalIndex = ref.getLocus().getStart() - ref.getWindow().getStart() + 1;
        final byte[] originalIndel = makeHaplotype(vc, refSeq, originalIndex, indelLength);

        // create a CIGAR string to represent the event
        ArrayList<CigarElement> elements = new ArrayList<CigarElement>();
        elements.add(new CigarElement(originalIndex, CigarOperator.M));
        elements.add(new CigarElement(indelLength, vc.isSimpleDeletion() ? CigarOperator.D : CigarOperator.I));
        elements.add(new CigarElement(refSeq.length - originalIndex, CigarOperator.M));
        Cigar originalCigar = new Cigar(elements);

        // left align the CIGAR
        Cigar newCigar = AlignmentUtils.leftAlignIndel(originalCigar, refSeq, originalIndel, 0, 0);

        // update if necessary and write
        if ( !newCigar.equals(originalCigar) && newCigar.numCigarElements() > 1 ) {
            int difference = originalIndex - newCigar.getCigarElement(0).getLength();
            VariantContext newVC = VariantContext.modifyLocation(vc, vc.getChr(), vc.getStart()-difference, vc.getEnd()-difference);
            //System.out.println("Moving record from " + vc.getChr()+":"+vc.getStart() + " to " + vc.getChr()+":"+(vc.getStart()-difference));

            int indelIndex = originalIndex-difference;
            byte[] newBases = new byte[indelLength];
            System.arraycopy((vc.isSimpleDeletion() ? refSeq : originalIndel), indelIndex, newBases, 0, indelLength);
            Allele newAllele = Allele.create(newBases, vc.isSimpleDeletion());
            newVC = updateAllele(newVC, newAllele, refSeq[indelIndex-1]);

            writer.add(newVC);
            return 1;
        } else {
            writer.add(vc);
            return 0;
        }
    }

    private static byte[] makeHaplotype(VariantContext vc, byte[] ref, int indexOfRef, int indelLength) {
        byte[] hap = new byte[ref.length + (indelLength * (vc.isSimpleDeletion() ? -1 : 1))];

        // add the bases before the indel
        System.arraycopy(ref, 0, hap, 0, indexOfRef);
        int currentPos = indexOfRef;

        // take care of the indel
        if ( vc.isSimpleDeletion() ) {
            indexOfRef += indelLength;
        } else {
            System.arraycopy(vc.getAlternateAllele(0).getBases(), 0, hap, currentPos, indelLength);
            currentPos += indelLength;
        }

        // add the bases after the indel
        System.arraycopy(ref, indexOfRef, hap, currentPos, ref.length - indexOfRef);

        return hap;
    }

    public static VariantContext updateAllele(VariantContext vc, Allele newAllele, Byte refBaseForIndel) {
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

        return new VariantContext(vc.getSource(), vc.getChr(), vc.getStart(), vc.getEnd(), alleleMap.values(), newGenotypes, vc.getNegLog10PError(), vc.filtersWereApplied() ? vc.getFilters() : null, vc.getAttributes(), refBaseForIndel);
    }
}
