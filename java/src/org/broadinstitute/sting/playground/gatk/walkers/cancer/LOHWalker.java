package org.broadinstitute.sting.playground.gatk.walkers.cancer;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.playground.utils.GenotypeLikelihoods;

import java.util.Map;
import java.util.HashMap;
import java.util.List;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: kcibul
 * Date: Jun 27, 2009
 * Time: 3:21:23 PM
 * To change this template use File | Settings | File Templates.
 */
public class LOHWalker extends LocusWalker<Integer, Integer> {
    @Argument(fullName = "tumor_sample_name", shortName = "s1", required = true, doc="Name of the tumor sample")
    public String tumorSampleName;

    @Argument(fullName = "normal_sample_name", shortName = "s2", required = true, doc="Name of the normal sample")
    public String normalSampleName;
                      
    public Integer reduceInit() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public Integer reduce(Integer value, Integer sum) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void initialize() {
        out.println("chrom\tposition\tnormal_btnb\tallele1\tallele2\tnormal_allele1_count\tnormal_allele2_count\tnormal_mar\ttumor_allele1_count\ttumor_allele2_count\ttumor_mar");

    }

    public static int MAX_INSERT_SIZE = 10000;

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        char upRef = Character.toUpperCase(ref);

        // TODO: should we be able to call mutations at bases where ref is N?
        if (upRef == 'N') { return -1; }

        List<SAMRecord> reads = context.getReads();

        LocusReadPile tumorReadPile = new LocusReadPile(upRef);
        LocusReadPile normalReadPile = new LocusReadPile(upRef);

        for ( int i = 0; i < reads.size(); i++ )
        {
            SAMRecord read = reads.get(i);

            if (read.getNotPrimaryAlignmentFlag() ||
                read.getDuplicateReadFlag() ||
                read.getReadUnmappedFlag() ||
                read.getMappingQuality() <= 0

//               || (read.getReadPairedFlag() && (!read.getProperPairFlag() || read.getInferredInsertSize() >= MAX_INSERT_SIZE))
                    ) {
                continue;
            }

            String rg = (String) read.getAttribute("RG");
            String sample = read.getHeader().getReadGroup(rg).getSample();

            int offset = context.getOffsets().get(i);

            char base = read.getReadString().charAt(offset);
            if (base == 'N' || base == 'n') { continue; }


            // TODO: build a pile of reads and offsets, then pass that into a
            // constructor for the normalGL class
            // that way, we can build a different pile of reads later on and extract the genotype
            if (normalSampleName.equals(sample)) {
                normalReadPile.add(read, offset);

            } else if (tumorSampleName.equals(sample)) {
                tumorReadPile.add(read, offset);
            } else {
                throw new RuntimeException("Unknown Sample Name: " + sample);
            }
        }

        normalReadPile.likelihoods.ApplyPrior(upRef,' ', -1); // apply std priors
        double normalBtnb = normalReadPile.likelihoods.LodVsNextBest();
        if (normalBtnb < 5.0) { return null; }

        String normalGt = normalReadPile.likelihoods.BestGenotype();
        char allele1 = normalGt.charAt(0);
        char allele2 = normalGt.charAt(1);
        if (allele2 == upRef) {
            char oldAllele1 = allele1;
            allele1 = allele2;
            allele2 = oldAllele1;
        }


        if (allele1 == upRef && allele1 == allele2) {
            return null;
        }
        
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("%s\t%d\t%1.2f\t%s\t%s\t",
            context.getContig(),
            context.getPosition(),
            normalBtnb,
            allele1,
            allele2));

        int normalAllele1Counts = normalReadPile.qualitySums.getCounts(allele1);
        int normalAllele2Counts = normalReadPile.qualitySums.getCounts(allele2);
        int tumorAllele1Counts = tumorReadPile.qualitySums.getCounts(allele1);
        int tumorAllele2Counts = tumorReadPile.qualitySums.getCounts(allele2);

        // allele ratio in the normal
        double normalAlleleRatio =
                (double) normalAllele2Counts / (double) (normalAllele1Counts+normalAllele2Counts) ;

        // allele ratio in the tumor
        double tumorAlleleRatio =
                (double) tumorAllele2Counts / (double) (tumorAllele1Counts+tumorAllele2Counts) ;


        sb.append(
                String.format("%d\t%d\t%1.2f\t%d\t%d\t%1.2f",
                                normalAllele1Counts,
                                normalAllele2Counts,
                                normalAlleleRatio,
                                tumorAllele1Counts,
                                tumorAllele2Counts,
                                tumorAlleleRatio));



//        // do the normal
//        for (char allele : new char[]{'A','C','G','T'}) {
//            QualitySums qs = normalReadPile.qualitySums;
//            sb.append(" ");
//            sb.append(qs.getCounts(allele)).append(" ");
////            sb.append(qs.getQualitySum(allele));
//        }
//
//        // do the tumor
//        for (char allele : new char[]{'A','C','G','T'}) {
//            QualitySums qs = tumorReadPile.qualitySums;
//            sb.append(" ");
//            sb.append(qs.getCounts(allele)).append(" ");
////            sb.append(qs.getQualitySum(allele));
//        }

        out.println(sb);
        return null;
    }
}
