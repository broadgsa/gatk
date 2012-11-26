package org.broadinstitute.sting.gatk.walkers.annotator;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.gatk.walkers.indels.PairHMMIndelErrorModel;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.*;

/**
 * The u-based z-approximation from the Mann-Whitney Rank Sum Test for the distance from the end of the read for reads with the alternate allele; if the alternate allele is only seen near the ends of reads this is indicative of error).
 * Note that the read position rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles.
 */
public class ReadPosRankSumTest extends RankSumTest implements StandardAnnotation {

    public List<String> getKeyNames() {
        return Arrays.asList("ReadPosRankSum");
    }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFInfoHeaderLine("ReadPosRankSum", 1, VCFHeaderLineType.Float, "Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias"));
    }

    protected void fillQualsFromPileup(final List<Allele> allAlleles,
                                       final int refLoc,
                                       final ReadBackedPileup pileup,
                                       final PerReadAlleleLikelihoodMap alleleLikelihoodMap,
                                       final List<Double> refQuals, final List<Double> altQuals) {

        if (alleleLikelihoodMap == null) {
            // use old UG SNP-based version if we don't have per-read allele likelihoods
            for ( final PileupElement p : pileup ) {
                if ( isUsableBase(p) ) {
                    int readPos = AlignmentUtils.calcAlignmentByteArrayOffset(p.getRead().getCigar(), p, 0, 0);

                    readPos = getFinalReadPosition(p.getRead(),readPos);

                    if ( allAlleles.get(0).equals(Allele.create(p.getBase(), true)) ) {
                        refQuals.add((double)readPos);
                    } else if ( allAlleles.contains(Allele.create(p.getBase()))) {
                        altQuals.add((double)readPos);
                    }
                }
            }
            return;
        }

        for (Map.Entry<GATKSAMRecord,Map<Allele,Double>> el : alleleLikelihoodMap.getLikelihoodReadMap().entrySet()) {
            final GATKSAMRecord read = el.getKey();
            final int offset = ReadUtils.getReadCoordinateForReferenceCoordinate( read.getSoftStart(), read.getCigar(), refLoc, ReadUtils.ClippingTail.RIGHT_TAIL, true );
            if ( offset == ReadUtils.CLIPPING_GOAL_NOT_REACHED )
                continue;
            int readPos = AlignmentUtils.calcAlignmentByteArrayOffset( read.getCigar(), offset, false, false, 0, 0 );
            final int numAlignedBases = AlignmentUtils.getNumAlignedBasesCountingSoftClips( read );
            if (readPos > numAlignedBases / 2)
                readPos = numAlignedBases - (readPos + 1);

//            int readPos = getOffsetFromClippedReadStart(el.getKey(), el.getKey().getOffset());
  //          readPos = getFinalReadPosition(el.getKey().getRead(),readPos);

            final Allele a = PerReadAlleleLikelihoodMap.getMostLikelyAllele(el.getValue());
            if (a.isNoCall())
                continue; // read is non-informative
            if (a.isReference())
                refQuals.add((double)readPos);
            else if (allAlleles.contains(a))
                altQuals.add((double)readPos);

        }
    }

    int getFinalReadPosition(GATKSAMRecord read, int initialReadPosition) {
        final int numAlignedBases = getNumAlignedBases(read);

        int readPos = initialReadPosition;
        if (initialReadPosition > numAlignedBases / 2) {
            readPos = numAlignedBases - (initialReadPosition + 1);
        }
        return readPos;

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
        for (int i = numStartClippedBases; i < unclippedReadBases.length; i++) {
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
        CigarElement last = c.getCigarElement(c.numCigarElements() - 1);

        int numEndClippedBases = 0;
        if (last.getOperator() == CigarOperator.H) {
            numEndClippedBases = last.getLength();
        }
        byte[] unclippedReadBases = read.getReadBases();
        byte[] unclippedReadQuals = read.getBaseQualities();

        // Do a stricter base clipping than provided by CIGAR string, since this one may be too conservative,
        // and may leave a string of Q2 bases still hanging off the reads.
        for (int i = unclippedReadBases.length - numEndClippedBases - 1; i >= 0; i--) {
            if (unclippedReadQuals[i] < PairHMMIndelErrorModel.BASE_QUAL_THRESHOLD)
                numEndClippedBases++;
            else
                break;
        }


        return numEndClippedBases;
    }

    int getOffsetFromClippedReadStart(SAMRecord read, int offset) {
        return offset - getNumClippedBasesAtStart(read);
    }
}
