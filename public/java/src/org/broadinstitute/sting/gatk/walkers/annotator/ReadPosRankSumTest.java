package org.broadinstitute.sting.gatk.walkers.annotator;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.walkers.genotyper.IndelGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.indels.PairHMMIndelErrorModel;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * The phred-scaled p-value (u-based z-approximation) from the Mann-Whitney Rank Sum Test for the distance from the end of the read for reads with the alternate allele; if the alternate allele is only seen near the ends of reads this is indicative of error).
 */
public class ReadPosRankSumTest extends RankSumTest {

    public List<String> getKeyNames() { return Arrays.asList("ReadPosRankSum"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("ReadPosRankSum", 1, VCFHeaderLineType.Float, "Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias")); }

    protected void fillQualsFromPileup(byte ref, byte alt, ReadBackedPileup pileup, List<Double> refQuals, List<Double> altQuals) {
        for ( final PileupElement p : pileup ) {
            if( isUsableBase(p) ) {
                int readPos = AlignmentUtils.calcAlignmentByteArrayOffset(p.getRead().getCigar(), p.getOffset(), 0, 0);
                final int numAlignedBases = AlignmentUtils.getNumAlignedBases(p.getRead());
                if( readPos > numAlignedBases / 2 ) {
                    readPos = numAlignedBases - ( readPos + 1 );
                }

                if ( p.getBase() == ref ) {
                    refQuals.add( (double)readPos );
                } else if ( p.getBase() == alt ) {
                    altQuals.add( (double)readPos );
                }
            }
        }
    }
    protected void fillIndelQualsFromPileup(ReadBackedPileup pileup, List<Double> refQuals, List<Double> altQuals) {
        // equivalent is whether indel likelihoods for reads corresponding to ref allele are more likely than reads corresponding to alt allele
        // to classify a pileup element as ref or alt, we look at the likelihood associated with the allele associated to this element.
        // A pileup element then has a list of pairs of form (Allele, likelihood of this allele).
        // To classify a pileup element as Ref or Alt, we look at the likelihood of corresponding alleles.
        // If likelihood of ref allele > highest likelihood of all alt alleles  + epsilon, then this pielup element is "ref"
        // otherwise  if highest alt allele likelihood is > ref likelihood + epsilon, then this pileup element it "alt"
        final HashMap<PileupElement,LinkedHashMap<Allele,Double>> indelLikelihoodMap = IndelGenotypeLikelihoodsCalculationModel.getIndelLikelihoodMap();
        for (final PileupElement p: pileup) {
            if (indelLikelihoodMap.containsKey(p)) {
                // retrieve likelihood information corresponding to this read
                LinkedHashMap<Allele,Double> el = indelLikelihoodMap.get(p);
                // by design, first element in LinkedHashMap was ref allele
                double refLikelihood=0.0, altLikelihood=Double.NEGATIVE_INFINITY;

                for (Allele a : el.keySet()) {

                    if (a.isReference())
                        refLikelihood =el.get(a);
                    else {
                        double like = el.get(a);
                        if (like >= altLikelihood)
                            altLikelihood = like;
                    }
                }

                int readPos = getOffsetFromClippedReadStart(p.getRead(), p.getOffset());
                final int numAlignedBases = getNumAlignedBases(p.getRead());

                int rp = readPos;
                if( readPos > numAlignedBases / 2 ) {
                    readPos = numAlignedBases - ( readPos + 1 );
                }
		//if (DEBUG) System.out.format("R:%s start:%d C:%s offset:%d rp:%d readPos:%d alignedB:%d\n",p.getRead().getReadName(),p.getRead().getAlignmentStart(),p.getRead().getCigarString(),p.getOffset(), rp, readPos, numAlignedBases);


                // if event is beyond span of read just return and don't consider this element. This can happen, for example, with reads
                // where soft clipping still left strings of low quality bases but these are later removed by indel-specific clipping.
               // if (readPos < -1)
                //    return;
                if (refLikelihood > (altLikelihood + INDEL_LIKELIHOOD_THRESH))   {
                    refQuals.add((double)readPos);
                    //if (DEBUG)  System.out.format("REF like: %4.1f, pos: %d\n",refLikelihood,readPos);
                }
                else if (altLikelihood > (refLikelihood + INDEL_LIKELIHOOD_THRESH))  {
                    altQuals.add((double)readPos);
		    //if (DEBUG)    System.out.format("ALT like: %4.1f, pos: %d\n",refLikelihood,readPos);

                }


            }
        }
    }

    int getNumClippedBasesAtStart(SAMRecord read) {
        // compute total number of clipped bases (soft or hard clipped)
        // check for hard clips (never consider these bases):
        final Cigar c = read.getCigar();
        final CigarElement first = c.getCigarElement(0);

        int numStartClippedBases = 0;
        if (first.getOperator() == CigarOperator.H) {
            numStartClippedBases = first.getLength();
        }
        byte[] unclippedReadBases = read.getReadBases();
        byte[] unclippedReadQuals = read.getBaseQualities();

        // Do a stricter base clipping than provided by CIGAR string, since this one may be too conservative,
        // and may leave a string of Q2 bases still hanging off the reads.
        for (int i=numStartClippedBases; i < unclippedReadBases.length; i++) {
            if (unclippedReadQuals[i] < PairHMMIndelErrorModel.BASE_QUAL_THRESHOLD)
                numStartClippedBases++;
            else
                break;

        }

        return numStartClippedBases;
    }

    int getNumAlignedBases(SAMRecord read) {
        return read.getReadLength() - getNumClippedBasesAtStart(read) - getNumClippedBasesAtEnd(read);
    }

    int getNumClippedBasesAtEnd(SAMRecord read) {
        // compute total number of clipped bases (soft or hard clipped)
        // check for hard clips (never consider these bases):
        final Cigar c = read.getCigar();
        CigarElement last = c.getCigarElement(c.numCigarElements()-1);

        int numEndClippedBases = 0;
        if (last.getOperator() == CigarOperator.H) {
            numEndClippedBases = last.getLength();
        }
        byte[] unclippedReadBases = read.getReadBases();
        byte[] unclippedReadQuals = read.getBaseQualities();

        // Do a stricter base clipping than provided by CIGAR string, since this one may be too conservative,
        // and may leave a string of Q2 bases still hanging off the reads.
        for (int i=unclippedReadBases.length-numEndClippedBases-1; i >= 0; i-- ){
            if (unclippedReadQuals[i] < PairHMMIndelErrorModel.BASE_QUAL_THRESHOLD)
                numEndClippedBases++;
            else
                break;
        }


        return numEndClippedBases;
    }

    int getOffsetFromClippedReadStart(SAMRecord read, int offset) {


        return offset -  getNumClippedBasesAtStart(read);
    }
}
