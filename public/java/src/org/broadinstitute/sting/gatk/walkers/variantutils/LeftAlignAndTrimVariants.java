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

package org.broadinstitute.sting.gatk.walkers.variantutils;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;

import java.util.*;

/**
 * Left-aligns indels from a variants file.
 *
 * <p>
 * LeftAlignAndTrimVariants is a tool that takes a VCF file and left-aligns the indels inside it.  The same indel can often be
 * placed at multiple positions and still represent the same haplotype.  While the standard convention with VCF is to
 * place an indel at the left-most position this doesn't always happen, so this tool can be used to left-align them.
 * Note that this tool cannot handle anything other than bi-allelic, simple indels.  Complex events are written out unchanged.
 * Optionally, the tool will also trim common bases from indels, leaving them with a minimum representation.
 *
 * <h3>Input</h3>
 * <p>
 * A variant set to left-align and trim.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A left-aligned VCF.
 * </p>
 *
 * <h3>Examples</h3>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T LeftAlignAndTrimVariants \
 *   --variant input.vcf \
 *   -o output.vcf
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
@Reference(window=@Window(start=-200,stop=200))    // WARNING: if this changes,MAX_INDEL_LENGTH needs to change as well!
public class LeftAlignAndTrimVariants extends RodWalker<Integer, Integer> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    /**
     * If this argument is set, bases common to all alleles will be removed, leaving only their minimal representation.
     */
    @Argument(fullName="trimAlleles", shortName="trim", doc="Trim alleles to remove bases common to all of them", required=false)
    protected boolean trimAlleles = false;

    /**
     * If this argument is set, split multiallelic records and left-align individual alleles.
     * If this argument is not set, multiallelic records are not attempted to left-align and will be copied as is.
     */
    @Argument(fullName="splitMultiallelics", shortName="split", doc="Split multiallelic records and left-align individual alleles", required=false)
    protected boolean splitMultiallelics = false;


    @Output(doc="File to which variants should be written")
    protected VariantContextWriter baseWriter = null;

    private VariantContextWriter writer;

    private static final int MAX_INDEL_LENGTH = 200; // needs to match reference window size!
    public void initialize() {
        String trackName = variantCollection.variants.getName();
        Set<String> samples = SampleUtils.getSampleListWithVCFHeader(getToolkit(), Arrays.asList(trackName));
        Map<String, VCFHeader> vcfHeaders = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), Arrays.asList(trackName));

        Set<VCFHeaderLine> headerLines = vcfHeaders.get(trackName).getMetaDataInSortedOrder();
        baseWriter.writeHeader(new VCFHeader(headerLines, samples));

        writer = VariantContextWriterFactory.sortOnTheFly(baseWriter, 200);
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        Collection<VariantContext> VCs = tracker.getValues(variantCollection.variants, context.getLocation());

        int changedSites = 0;
        for ( final VariantContext vc : VCs ) {
            // split first into biallelics, and optionally trim alleles to minimal representation
            Pair<VariantContext,Integer> result = new Pair<VariantContext, Integer>(vc,0); // default value
            if (splitMultiallelics) {
                final List<VariantContext> vcList = GATKVariantContextUtils.splitVariantContextToBiallelics( vc);
                for (final VariantContext biallelicVC: vcList) {
                    final VariantContext v = (trimAlleles ? GATKVariantContextUtils.trimAlleles(biallelicVC,true,true):biallelicVC);
                    result = alignAndWrite(v, ref);
                    writer.add(result.first);
                    changedSites += result.second;
                }
            }
            else {
                if (trimAlleles)
                    result = alignAndWrite(GATKVariantContextUtils.trimAlleles(vc,true,true), ref);
                else
                    result = alignAndWrite(vc,ref);
                writer.add(result.first);
                changedSites += result.second;

            }

        }

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

    /**
     * Main routine workhorse. By definitio, it will only take biallelic vc's. Splitting into multiple alleles has to be
     * handled by calling routine.
     * @param vc                  Input VC with variants to left align
     * @param ref                 Reference context
     * @return                    # of records left-aligned (0 or 1) and new VC.
     */
    @Requires({"vc != null","ref != null", "vc.isBiallelic() == true","ref.getBases().length>=2*MAX_INDEL_LENGTH+1"})
    @Ensures({"result != null","result.first != null", "result.second >=0"})
    protected static Pair<VariantContext,Integer>  alignAndWrite(final VariantContext vc, final ReferenceContext ref) {

        final Pair<VariantContext, Integer> retValue =  new Pair<VariantContext, Integer>(vc,0);
        if (!vc.isIndel() || vc.isComplexIndel() ) {
            return retValue;
        }

        // get the indel length
        final int indelLength;
        if ( vc.isSimpleDeletion() )
            indelLength = vc.getReference().length() - 1;
        else
            indelLength = vc.getAlternateAllele(0).length() - 1;

        if ( indelLength > MAX_INDEL_LENGTH )
            return retValue;

         if (vc.getReference().getBases()[0] != vc.getAlternateAllele(0).getBases()[0])
            return retValue;

        final byte[] refSeq = ref.getBases();

        // create an indel haplotype.
        //
        final int originalIndex = vc.getStart() - ref.getWindow().getStart() + 1;
        if (originalIndex < 0 || originalIndex >= ref.getBases().length)
            return retValue;

        final byte[] originalIndel = makeHaplotype(vc, refSeq, originalIndex, indelLength);

        // create a CIGAR string to represent the event
        ArrayList<CigarElement> elements = new ArrayList<CigarElement>();
        elements.add(new CigarElement(originalIndex, CigarOperator.M));
        elements.add(new CigarElement(indelLength, vc.isSimpleDeletion() ? CigarOperator.D : CigarOperator.I));
        elements.add(new CigarElement(refSeq.length - originalIndex, CigarOperator.M));
        Cigar originalCigar = new Cigar(elements);

        // left align the CIGAR
        Cigar newCigar = AlignmentUtils.leftAlignIndel(originalCigar, refSeq, originalIndel, 0, 0, true);

        // update if necessary and write
        if ( !newCigar.equals(originalCigar) && newCigar.numCigarElements() > 1 ) {
            int difference = originalIndex - newCigar.getCigarElement(0).getLength();
            VariantContext newVC = new VariantContextBuilder(vc).start(vc.getStart()-difference).stop(vc.getEnd()-difference).make();
            //System.out.println("Moving record from " + vc.getChr()+":"+vc.getStart() + " to " + vc.getChr()+":"+(vc.getStart()-difference));

            final int indelIndex = originalIndex-difference;
            final byte[] newBases = new byte[indelLength + 1];
            newBases[0] = refSeq[indelIndex-1];
            System.arraycopy((vc.isSimpleDeletion() ? refSeq : originalIndel), indelIndex, newBases, 1, indelLength);
            final Allele newAllele = Allele.create(newBases, vc.isSimpleDeletion());
            newVC = updateAllele(newVC, newAllele);
            // overwrite default return value with new left-aligned VC
            retValue.first = newVC;
            retValue.second = 1;

        }
        return retValue;
    }

    /**
     * Make a haplotype from a given alt allele, using bases in input reference, index of an input reference
     * @param vc                                Input VC - will use only alt allele from it
     * @param ref                               Ref bases
     * @param indexOfRef                        Index in ref where to create indel
     * @param indelLength                       Indel length
     * @return
     */
    @Requires({"vc != null","ref != null", "indexOfRef +indelLength < ref.length", "vc.getNAlleles() == 2"})
    @Ensures("result != null")
    private static byte[] makeHaplotype(VariantContext vc, byte[] ref, int indexOfRef, int indelLength) {
        byte[] hap = new byte[ref.length + (indelLength * (vc.isSimpleDeletion() ? -1 : 1))];

        // add the bases before the indel
        System.arraycopy(ref, 0, hap, 0, indexOfRef);
        int currentPos = indexOfRef;

        // take care of the indel
        if ( vc.isSimpleDeletion() ) {
            indexOfRef += indelLength;
        } else {
            System.arraycopy(vc.getAlternateAllele(0).getBases(), 1, hap, currentPos, indelLength);
            currentPos += indelLength;
        }

        // add the bases after the indel
        System.arraycopy(ref, indexOfRef, hap, currentPos, ref.length - indexOfRef);

        return hap;
    }

    public static VariantContext updateAllele(final VariantContext vc, final Allele newAllele) {
        // create a mapping from original allele to new allele
        HashMap<Allele, Allele> alleleMap = new HashMap<Allele, Allele>(vc.getAlleles().size());
        if ( newAllele.isReference() ) {
            alleleMap.put(vc.getReference(), newAllele);
            alleleMap.put(vc.getAlternateAllele(0), Allele.create(newAllele.getBases()[0], false));
        } else {
            alleleMap.put(vc.getReference(), Allele.create(newAllele.getBases()[0], true));
            alleleMap.put(vc.getAlternateAllele(0), newAllele);
        }

        // create new Genotype objects
        GenotypesContext newGenotypes = GenotypesContext.create(vc.getNSamples());
        for ( final Genotype genotype : vc.getGenotypes() ) {
            List<Allele> newAlleles = new ArrayList<Allele>();
            for ( Allele allele : genotype.getAlleles() ) {
                Allele newA = alleleMap.get(allele);
                if ( newA == null )
                    newA = Allele.NO_CALL;
                newAlleles.add(newA);
            }
            newGenotypes.add(new GenotypeBuilder(genotype).alleles(newAlleles).make());
        }

        return new VariantContextBuilder(vc).alleles(alleleMap.values()).genotypes(newGenotypes).make();
    }
}
