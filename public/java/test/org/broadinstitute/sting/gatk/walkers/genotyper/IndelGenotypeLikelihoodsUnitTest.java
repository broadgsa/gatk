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
    
    final int contigStart = 1;
    final int contigStop = 10;
    final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, contigStop-contigStart+1);
    final GATKSAMReadGroupRecord artificialGATKRG = new GATKSAMReadGroupRecord("synthetic");
    final String artificialContig = "chr1";
    final int artificialContigIndex = 0;
    final String artificialReadName = "synth";
    final int artificialRefStart = 1;
    final int artificialMappingQuality = 60;
    Map<String, SAMReadGroupRecord> sample2RG = new HashMap<String, SAMReadGroupRecord>();

    final String refBases = "AGGATACTGT";
    final String SAMPLE_PREFIX = "sample";
    
    List<String> sampleNames = new ArrayList<String>();
    final int nSamples = 1;
    final int numReadsPerAllele = 10;

    List<SAMReadGroupRecord> sampleRGs;

    private String sampleName(int i) { return sampleNames.get(i); }
    private SAMReadGroupRecord sampleRG(String name) { return sample2RG.get(name); }


    final UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();
    final Logger logger = Logger.getLogger(Walker.class);
    final GenomeLocParser genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    final IndelGenotypeLikelihoodsCalculationModel model = new IndelGenotypeLikelihoodsCalculationModel(UAC,logger);
    final int offset = 5;
    final GenomeLoc loc = genomeLocParser.createGenomeLoc(artificialContig,offset,offset);
    final GenomeLoc window = genomeLocParser.createGenomeLoc(artificialContig,artificialRefStart,10);
    final ReferenceContext referenceContext = new ReferenceContext(genomeLocParser,loc,window,this.refBases.getBytes());

    @BeforeSuite
    public void before() {
        sampleRGs = new ArrayList<SAMReadGroupRecord>();

        for ( int i = 0; i < nSamples; i++ ) {
            sampleNames.add(String.format("%s%04d", SAMPLE_PREFIX, i));
            SAMReadGroupRecord rg = createRG(sampleName(i));
            sampleRGs.add(rg);
            sample2RG.put(sampleName(i), rg);
        }

    }
    @Test
    public void testBasicConsensusCounts() {
        // 4 inserted bases, min cnt = 10
        String altBases = "CCTCCTGAGA";
        int eventLength = 4;
        List<Allele> alleles = getConsensusAlleles(eventLength,true,10,0.1, altBases);
        
        Assert.assertEquals(alleles.size(),2);
        Assert.assertEquals(alleles.get(1).getBaseString(), altBases.substring(0,eventLength));



        //altBases = "CCTCMTGAGA";

        eventLength = 3;
        alleles = getConsensusAlleles(eventLength,false,10,0.1, altBases);
        Assert.assertEquals(alleles.size(),2);
        Assert.assertEquals(alleles.get(0).getBaseString(), refBases.substring(offset,offset+eventLength));

        // same with min Reads = 11
        alleles = getConsensusAlleles(eventLength,false,11,0.1, altBases);
        Assert.assertEquals(alleles.size(),0);

        // increase required fraction per sample to just below threshold
        alleles = getConsensusAlleles(eventLength,false,10,0.49999, altBases);
        Assert.assertEquals(alleles.size(),2);
        alleles = getConsensusAlleles(eventLength,false,10,0.5001, altBases);
        Assert.assertEquals(alleles.size(),0);
    }
    
    private List<Allele> getConsensusAlleles(int eventLength, boolean isInsertion, int minCnt, double minFraction, String altBases) {
        final ConsensusAlleleCounter counter = new ConsensusAlleleCounter(genomeLocParser, true, minCnt, minFraction);
        return counter.computeConsensusAlleles(referenceContext,getContextFromAlleles(eventLength, isInsertion, altBases), AlignmentContextUtils.ReadOrientation.COMPLETE);

    }
    private Map<String,AlignmentContext> getContextFromAlleles(int eventLength, boolean isInsertion, String altBases) {
    //    RefMetaDataTracker tracker = new RefMetaDataTracker(null,referenceContext);

        
        ArrayList vcAlleles = new ArrayList<Allele>();
        Allele refAllele, altAllele;
        if (isInsertion) {
            refAllele = Allele.create(Allele.NULL_ALLELE_STRING, true);
            altAllele = Allele.create(altBases.substring(0,eventLength), false);
        }
        else {
            refAllele =Allele.create(refBases.substring(offset,offset+eventLength),true);
            altAllele = Allele.create(Allele.NULL_ALLELE_STRING, false);
        }
   
        int stop = loc.getStart();
        vcAlleles.add(refAllele);
        vcAlleles.add(altAllele);

        final VariantContextBuilder builder = new VariantContextBuilder().source("");
        builder.loc(loc.getContig(), loc.getStart(), stop);
        builder.alleles(vcAlleles);
        builder.referenceBaseForIndel(referenceContext.getBase());
        builder.noGenotypes();
        
        final VariantContext vc = builder.make();

        Map<String,AlignmentContext> contexts = new HashMap<String,AlignmentContext>();

        for (String sample: sampleNames) {
            AlignmentContext context = new AlignmentContext(loc, generateRBPForVariant(loc,vc, altBases, numReadsPerAllele, sample));
            contexts.put(sample,context);
            
        }

        return contexts;
    }

    private SAMReadGroupRecord createRG(String name) {
        SAMReadGroupRecord rg = new SAMReadGroupRecord(name);
        rg.setPlatform("ILLUMINA");
        rg.setSample(name);
        return rg;
    }
    private ReadBackedPileup generateRBPForVariant( GenomeLoc loc, VariantContext vc, String altBases,
                                                    int numReads, String sample) {
        List<PileupElement> pileupElements = new ArrayList<PileupElement>();
        int readStart = contigStart;
        int offset = (contigStop-contigStart+1)/2;
        int refAlleleLength = 0;
        int readCounter = 0;
        for (Allele allele: vc.getAlleles()) {
            if (allele.isReference())
                refAlleleLength = allele.getBases().length;

            int alleleLength = allele.getBases().length;

            for ( int d = 0; d < numReads; d++ ) {
                byte[] readBases = trueHaplotype(allele, offset, refAlleleLength);
                byte[] readQuals = new byte[readBases.length];
                Arrays.fill(readQuals,(byte)50);

                GATKSAMRecord read = new GATKSAMRecord(header);
                read.setBaseQualities(readQuals);
                read.setReadBases(readBases);
                read.setReadName(artificialReadName+readCounter++);

                boolean isBeforeDeletion = false, isBeforeInsertion = false;
                if (allele.isReference())
                    read.setCigarString(readBases.length + "M");
                else {
                    isBeforeDeletion = alleleLength<refAlleleLength;
                    isBeforeInsertion = alleleLength>refAlleleLength;
                    read.setCigarString(offset+"M"+ alleleLength + (isBeforeDeletion?"D":"I") +
                            (readBases.length-offset)+"M");
                }

                int eventLength = (isBeforeDeletion?refAlleleLength:(isBeforeInsertion?alleleLength:0));
                read.setReadPairedFlag(false);
                read.setAlignmentStart(readStart);
                read.setMappingQuality(artificialMappingQuality);
                read.setReferenceName(loc.getContig());
                read.setReadNegativeStrandFlag(false);
                read.setAttribute("RG", sampleRG(sample).getReadGroupId());


                pileupElements.add(new PileupElement(read,offset,false,isBeforeDeletion, false, isBeforeInsertion,false,false,altBases.substring(0,alleleLength),eventLength));
            }
        }

        return new ReadBackedPileupImpl(loc,pileupElements);
    }

    byte[] trueHaplotype(Allele allele, int offset, int refAlleleLength) {
        // create haplotype based on a particular allele
        String prefix = refBases.substring(offset);
        String alleleBases = new String(allele.getBases());
        String postfix = refBases.substring(offset+refAlleleLength,refBases.length());
        
        return (prefix+alleleBases+postfix).getBytes();



    }
}
