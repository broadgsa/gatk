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
package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

import java.util.*;


public class ArtificialReadPileupTestProvider {
    final int contigStart = 1;
    final int contigStop = 10;
    final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, contigStop - contigStart + 1);
//    final GATKSAMReadGroupRecord artificialGATKRG = new GATKSAMReadGroupRecord("synthetic");
    final String artificialContig = "chr1";
  //  final int artificialContigIndex = 0;
    final String artificialReadName = "synth";
    final int artificialRefStart = 1;
    final int artificialMappingQuality = 60;
    Map<String, SAMReadGroupRecord> sample2RG = new HashMap<String, SAMReadGroupRecord>();
    List<SAMReadGroupRecord> sampleRGs;

    final String refBases = "AGGATACTGT";
    List<String> sampleNames = new ArrayList<String>();
    private String sampleName(int i) { return sampleNames.get(i); }
    private SAMReadGroupRecord sampleRG(String name) { return sample2RG.get(name); }
    public final int offset = 5;
    public final GenomeLocParser genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    public final GenomeLoc loc = genomeLocParser.createGenomeLoc(artificialContig,offset,offset);
    public final GenomeLoc window = genomeLocParser.createGenomeLoc(artificialContig,artificialRefStart,10);
    public final ReferenceContext referenceContext = new ReferenceContext(genomeLocParser,loc,window,this.refBases.getBytes());

    byte BASE_QUAL = 50;

    public ArtificialReadPileupTestProvider(final int numSamples, final String SAMPLE_PREFIX) {
        sampleRGs = new ArrayList<SAMReadGroupRecord>();

        for ( int i = 0; i < numSamples; i++ ) {
            sampleNames.add(String.format("%s%04d", SAMPLE_PREFIX, i));
            SAMReadGroupRecord rg = createRG(sampleName(i));
            sampleRGs.add(rg);
            sample2RG.put(sampleName(i), rg);
        }

    }

    public ArtificialReadPileupTestProvider(final int numSamples, final String SAMPLE_PREFIX, final byte q) {
        this(numSamples,SAMPLE_PREFIX);
        BASE_QUAL = q;
    }
    public List<String> getSampleNames() {
        return sampleNames;
    }
    public byte getRefByte() {
        return referenceContext.getBase();
    }

    public ReferenceContext getReferenceContext()   { return referenceContext;}
    public GenomeLocParser getGenomeLocParser()     { return genomeLocParser; }

    public Map<String,AlignmentContext> getAlignmentContextFromAlleles(int eventLength, String altBases, int[] numReadsPerAllele) {
        return getAlignmentContextFromAlleles(eventLength, altBases, numReadsPerAllele, false, BASE_QUAL);
    }
    public Map<String,AlignmentContext> getAlignmentContextFromAlleles(int eventLength, String altBases, int[] numReadsPerAllele,
                                                                       boolean addBaseErrors, int phredScaledBaseErrorRate) {
        //    RefMetaDataTracker tracker = new RefMetaDataTracker(null,referenceContext);

        ArrayList<Allele> vcAlleles = new ArrayList<Allele>();
        String refBase = refBases.substring(offset,offset+1);    // referenceContext.getBase()?
        Allele refAllele, altAllele;
        if (eventLength == 0)  {
            // SNP case
            refAllele = Allele.create(refBase,true);
            altAllele = Allele.create(altBases.substring(0,1), false);

        } else if (eventLength>0){
            // insertion
            refAllele = Allele.create(refBase,true);
            altAllele = Allele.create(refBase + altBases.substring(0,eventLength), false);
        }
        else {
            // deletion
            refAllele = Allele.create(refBases.substring(offset,offset+Math.abs(eventLength)),true);
            altAllele = Allele.create(refBase, false);
        }
        int stop = loc.getStart();
        vcAlleles.add(refAllele);
        vcAlleles.add(altAllele);

        final VariantContextBuilder builder = new VariantContextBuilder().source("");
        builder.loc(loc.getContig(), loc.getStart(), stop);
        builder.alleles(vcAlleles);
        builder.noGenotypes();

        final VariantContext vc = builder.make();

        Map<String,AlignmentContext> contexts = new HashMap<String,AlignmentContext>();

        for (String sample: sampleNames) {
            AlignmentContext context = new AlignmentContext(loc, generateRBPForVariant(loc,vc, altBases, numReadsPerAllele, sample, addBaseErrors, phredScaledBaseErrorRate));
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
                                                    int[] numReadsPerAllele, String sample, boolean addErrors, int phredScaledErrorRate) {
        List<PileupElement> pileupElements = new ArrayList<PileupElement>();
        int readStart = contigStart;
        int offset = (contigStop-contigStart+1)/2;
        int refAlleleLength = 0;
        int readCounter = 0;
        int alleleCounter = 0;
        for (Allele allele: vc.getAlleles()) {
            if (allele.isReference())
                refAlleleLength = allele.getBases().length;

            int alleleLength = allele.getBases().length;

            for ( int d = 0; d < numReadsPerAllele[alleleCounter]; d++ ) {
                byte[] readBases = trueHaplotype(allele, offset, refAlleleLength);
                if (addErrors)
                    addBaseErrors(readBases, phredScaledErrorRate);

                byte[] readQuals = new byte[readBases.length];
                Arrays.fill(readQuals, (byte)phredScaledErrorRate);

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
                    if (isBeforeDeletion || isBeforeInsertion)
                        read.setCigarString(offset+"M"+ alleleLength + (isBeforeDeletion?"D":"I") +
                            (readBases.length-offset)+"M");
                    else // SNP case
                        read.setCigarString(readBases.length+"M");
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
            alleleCounter++;
        }

        return new ReadBackedPileupImpl(loc,pileupElements);
    }

    private byte[] trueHaplotype(Allele allele, int offset, int refAlleleLength) {
        // create haplotype based on a particular allele
        String prefix = refBases.substring(offset);
        String alleleBases = new String(allele.getBases());
        String postfix = refBases.substring(offset+refAlleleLength,refBases.length());

        return (prefix+alleleBases+postfix).getBytes();



    }

    private void addBaseErrors(final byte[] readBases, final int phredScaledErrorRate) {
        double errorProbability = QualityUtils.qualToErrorProb((byte)phredScaledErrorRate);

        for (int k=0; k < readBases.length; k++) {
            if (GenomeAnalysisEngine.getRandomGenerator().nextDouble() < errorProbability) {
                // random offset
                int offset = BaseUtils.simpleBaseToBaseIndex(readBases[k]);          //0..3
                offset += (GenomeAnalysisEngine.getRandomGenerator().nextInt(3)+1);  // adds 1,2 or 3
                offset %= 4;
                readBases[k] = BaseUtils.baseIndexToSimpleBase(offset);

            }

        }

    }
}
