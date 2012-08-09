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
    final String refBases = "ACAGAGCTGACCCTCCCTCCCCTCTCCCAGTGCAACAGCACGGGCGGCGACTGCTTTTACCGAGGCTACACGTCAGGCGTGGCGGCTGTCCAGGACTGGTACCACTTCCACTATGTGGATCTCTGCTGAGGACCAGGAAAGCCAGCACCCGCAGAGACTCTTCCCCAGTGCTCCATACGATCACCATTCTCTGCAGAAGGTCAGACGTCACTGGTGGCCCCCCAGCCTCCTCAGCAGGGAAGGATACTGTCCCGCAGATGAGATGAGCGAGAGCCGCCAGACCCACGTGACGCTGCACGACATCGACCCTCAGGCCTTGGACCAGCTGGTGCAGTTTGCCTACACGGCTGAGATTGTGGTGGGCGAGGGC";
    final int contigStart = 1;
    final int contigStop = refBases.length();
    final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, contigStop - contigStart + 1);
//    final GATKSAMReadGroupRecord artificialGATKRG = new GATKSAMReadGroupRecord("synthetic");
    final String artificialContig = "chr1";
  //  final int artificialContigIndex = 0;
    final String artificialReadName = "synth";
    final int artificialRefStart = 1;
    final int artificialMappingQuality = 60;
    Map<String, SAMReadGroupRecord> sample2RG = new HashMap<String, SAMReadGroupRecord>();
    List<SAMReadGroupRecord> sampleRGs;
    List<String> sampleNames = new ArrayList<String>();
    private String sampleName(int i) { return sampleNames.get(i); }
    private SAMReadGroupRecord sampleRG(String name) { return sample2RG.get(name); }
    public final int locStart = 105; // start position where we desire artificial variant
    private final int readLength = 10; // desired read length in pileup
    public final int readOffset = 4;
    private final int readStart = locStart - readOffset;
    public final GenomeLocParser genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    public final GenomeLoc loc = genomeLocParser.createGenomeLoc(artificialContig,locStart,locStart);
    public final GenomeLoc window = genomeLocParser.createGenomeLoc(artificialContig,locStart-100,locStart+100);
    public final String windowBases = refBases.substring(locStart-100-1,locStart+100);
    public final ReferenceContext referenceContext = new ReferenceContext(genomeLocParser,loc,window,windowBases.getBytes());

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
    public Map<String,AlignmentContext> getAlignmentContextFromAlleles(final int eventLength,
                                                                       final String altBases,
                                                                       final int[] numReadsPerAllele,
                                                                       final boolean addBaseErrors,
                                                                       final int phredScaledBaseErrorRate) {
        final String refChar = new String(new byte[]{referenceContext.getBase()});

        String refAllele, altAllele;
        if (eventLength == 0)  {
            // SNP case
            refAllele = refChar;
            altAllele = altBases.substring(0,1);

        } else if (eventLength>0){
            // insertion
            refAllele = refChar;
            altAllele = refChar+altBases/*.substring(0,eventLength)*/;
        }
        else {
            // deletion
            refAllele = new String(referenceContext.getForwardBases()).substring(0,Math.abs(eventLength)+1);
            altAllele = refChar;
        }

        Map<String,AlignmentContext> contexts = new HashMap<String,AlignmentContext>();

        for (String sample: sampleNames) {
            AlignmentContext context = new AlignmentContext(loc, generateRBPForVariant(loc, refAllele, altAllele, altBases, numReadsPerAllele, sample, addBaseErrors, phredScaledBaseErrorRate));
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

    private ReadBackedPileup generateRBPForVariant( GenomeLoc loc, String refAllele, String altAllele, String altBases,
                                                    int[] numReadsPerAllele, String sample, boolean addErrors, int phredScaledErrorRate) {
        List<PileupElement> pileupElements = new ArrayList<PileupElement>();
        final int refAlleleLength = refAllele.length();

        pileupElements.addAll(createPileupElements(refAllele, loc, numReadsPerAllele[0], sample, readStart, altBases, addErrors, phredScaledErrorRate, refAlleleLength, true));
        pileupElements.addAll(createPileupElements(altAllele, loc, numReadsPerAllele[1], sample, readStart, altBases, addErrors, phredScaledErrorRate, refAlleleLength, false));
        return new ReadBackedPileupImpl(loc,pileupElements);
    }

    private List<PileupElement> createPileupElements(String allele, GenomeLoc loc, int numReadsPerAllele, String sample, int readStart, String altBases, boolean addErrors, int phredScaledErrorRate, int refAlleleLength, boolean isReference) {

        int alleleLength = allele.length();
        List<PileupElement> pileupElements = new ArrayList<PileupElement>();

        int readCounter = 0;
        for ( int d = 0; d < numReadsPerAllele; d++ ) {
            byte[] readBases = trueHaplotype(allele, refAlleleLength, readLength);
            if (addErrors)
                addBaseErrors(readBases, phredScaledErrorRate);

            byte[] readQuals = new byte[readBases.length];
            Arrays.fill(readQuals, (byte)phredScaledErrorRate);

            GATKSAMRecord read = new GATKSAMRecord(header);
            read.setBaseQualities(readQuals);
            read.setReadBases(readBases);
            read.setReadName(artificialReadName+readCounter++);

            boolean isBeforeDeletion = alleleLength<refAlleleLength;
            boolean isBeforeInsertion = alleleLength>refAlleleLength;

            int eventLength = alleleLength - refAlleleLength;
            if (isReference)
                read.setCigarString(readBases.length + "M");
            else {
                if (isBeforeDeletion || isBeforeInsertion)
                    read.setCigarString((readOffset+1)+"M"+ Math.abs(eventLength) + (isBeforeDeletion?"D":"I") +
                            (readBases.length-readOffset)+"M");
                else // SNP case
                    read.setCigarString(readBases.length+"M");
            }

            read.setReadPairedFlag(false);
            read.setAlignmentStart(readStart);
            read.setMappingQuality(artificialMappingQuality);
            read.setReferenceName(loc.getContig());
            read.setReadNegativeStrandFlag(false);
            read.setAttribute("RG", sampleRG(sample).getReadGroupId());


            pileupElements.add(new PileupElement(read,readOffset,false,isBeforeDeletion, false, isBeforeInsertion,false,false,altBases,Math.abs(eventLength)));
        }

        return pileupElements;
    }

    /**
     * Create haplotype with desired allele and reference context
     * @param allele                             Desired allele string
     * @param refAlleleLength                    Length of reference allele.
     * @param desiredLength                      Desired haplotype length
     * @return                                   String with haplotype formed by (prefix)+allele bases + postfix
     */
    private byte[] trueHaplotype(final String allele, final int refAlleleLength, final int desiredLength) {
        // create haplotype based on a particular allele
        final int startIdx= locStart - readOffset-1;

        final String prefix = refBases.substring(startIdx, locStart-1);
        final String postfix = refBases.substring(locStart+refAlleleLength-1,startIdx + desiredLength);

        return (prefix+allele+postfix).getBytes();
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
