package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.genotyper.IndelGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

import java.util.*;
import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;

import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: 3/22/12
 * Time: 11:24 AM
 * To change this template use File | Settings | File Templates.
 */
public class IndelGenotypeLikelihoodsUnitTest extends BaseTest {
    
    final int nSamples = 1;
    final int[] numReadsPerAllele = new int[]{10,10};
    final String SAMPLE_PREFIX = "sample";


    final UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();
    final Logger logger = Logger.getLogger(Walker.class);
    final IndelGenotypeLikelihoodsCalculationModel model = new IndelGenotypeLikelihoodsCalculationModel(UAC,logger);

    ArtificialReadPileupTestProvider pileupProvider;

   @BeforeSuite
    public void before() {
        pileupProvider = new ArtificialReadPileupTestProvider(nSamples, SAMPLE_PREFIX);
    }

    @Test
    public void testBasicConsensusCounts() {
        // 4 inserted bases, min cnt = 10
        String altBases = "CCTC";
        int eventLength = 4;
        List<Allele> alleles = getConsensusAlleles(eventLength,true,10,0.1, altBases);
        
        Assert.assertEquals(alleles.size(),2);
        Assert.assertEquals(alleles.get(1).getBaseString().substring(1), altBases.substring(0,eventLength));


        // test deletions
        eventLength = 3;
        alleles = getConsensusAlleles(eventLength,false,10,0.1, altBases);
        Assert.assertEquals(alleles.size(),2);
        Assert.assertEquals(alleles.get(0).getBaseString().substring(1,eventLength), new String(pileupProvider.getReferenceContext().getForwardBases()).substring(1,eventLength));

        // same with min Reads = 11
        alleles = getConsensusAlleles(eventLength,false,11,0.1, altBases);
        Assert.assertEquals(alleles.size(),0);

        // increase required fraction per sample to just below threshold
        alleles = getConsensusAlleles(eventLength,false,10,0.49999, altBases);
        Assert.assertEquals(alleles.size(),2);
        alleles = getConsensusAlleles(eventLength,false,10,0.5001, altBases);
        Assert.assertEquals(alleles.size(),0);

        // test N's in insertions
        altBases = "CCTC";
        eventLength = 4;
        alleles = getConsensusAlleles(eventLength,true,10,0.1, altBases);

        Assert.assertEquals(alleles.size(),2);
        Assert.assertEquals(alleles.get(1).getBaseString().substring(1,eventLength+1), altBases);

        altBases = "CCTCN";
        eventLength = 5;
        alleles = getConsensusAlleles(eventLength,true,10,0.1, altBases);

        Assert.assertEquals(alleles.size(),0);

    }
    
    private List<Allele> getConsensusAlleles(int eventLength, boolean isInsertion, int minCnt, double minFraction, String altBases) {
        final ConsensusAlleleCounter counter = new ConsensusAlleleCounter(pileupProvider.genomeLocParser, true, minCnt, minFraction);
        return counter.computeConsensusAlleles(pileupProvider.referenceContext,
                pileupProvider.getAlignmentContextFromAlleles(isInsertion?eventLength:-eventLength,altBases,numReadsPerAllele),
                AlignmentContextUtils.ReadOrientation.COMPLETE);

    }
}
