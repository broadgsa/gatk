package org.broadinstitute.sting.playground.fourbasecaller;

import java.io.File;
import java.util.Vector;
import java.lang.Math;

import org.broadinstitute.sting.playground.illumina.FirecrestFileParser;
import org.broadinstitute.sting.playground.illumina.FourIntensity;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleFactory1D;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;

public class FourBaseCaller {
    public static void main(String[] argv) {
        // Parse args
        File FIRECREST_DIR = new File(argv[0]);
        int LANE = Integer.valueOf(argv[1]);
        File SAM_OUT = new File(argv[2]);
        int CYCLE_START = Integer.valueOf(argv[3]);
        int CYCLE_STOP = Integer.valueOf(argv[4]);

        String IC_REF_PATH = "/seq/references/Synthetic_internal_controls_set1/v0/Synthetic_internal_controls_set1.fasta";

        // Load ic reference in an ugly way
        ManyNucleotideSequences fw_controls = new ManyNucleotideSequences(IC_REF_PATH);
        ManyNucleotideSequences rc_controls = new ManyNucleotideSequences(IC_REF_PATH);
        rc_controls.reverseComplement();
        ManyNucleotideSequences controls = new ManyNucleotideSequences();
        controls.addAll(fw_controls);
        controls.addAll(rc_controls);

        // Loop through firecrest data, find ICs, and compute signal means
        Vector<FourIntensity[]> recoveredICs = new Vector<FourIntensity[]>();
        Vector<Integer> icids = new Vector<Integer>();

        FirecrestFileParser ffp = new FirecrestFileParser(FIRECREST_DIR, LANE, CYCLE_START, CYCLE_STOP);

        NucleotideChannelMeans[] cmeans = new NucleotideChannelMeans[controls.get(0).size()];
        for (int i = 0; i < cmeans.length; i++) {
            cmeans[i] = new NucleotideChannelMeans();
        }

        int numReads = 0;
        while (ffp.hasNext()) {
            if (numReads % 1000000 == 0) {
                System.err.println("Processed " + numReads + " reads for parameters.");
            }
            FourIntensity[] intensities = ffp.next().getIntensities();
            NucleotideSequence rawread = toNucleotideSequence(intensities);

            int[] editDistances = new int[controls.size()];

            for (int icTemplateIndex = 0; icTemplateIndex < controls.size(); icTemplateIndex++) {
                editDistances[icTemplateIndex] = rawread.editDistance(controls.get(icTemplateIndex), 36);

                if (editDistances[icTemplateIndex] <= 8) {
                    recoveredICs.add(intensities);
                    icids.add(icTemplateIndex);

                    for (int cycle = 0; cycle < ((controls.get(0).size() < intensities.length) ? controls.get(0).size() : intensities.length); cycle++) {
                        Nucleotide nuc = controls.get(icTemplateIndex).get(cycle);
                        FourIntensity sig = intensities[cycle];

                        cmeans[cycle].add(nuc, sig);
                    }

                    System.err.print("Found " + icids.size() + " ICs in " + numReads + " reads " + (((double) icids.size())/((double) numReads)) + ".\r");

                    break;
                }
            }

            numReads++;
        }

        // Go through the data again and compute signal covariances
        NucleotideChannelCovariances[] ccov = new NucleotideChannelCovariances[controls.get(0).size()];
        for (int i = 0; i < cmeans.length; i++) {
            ccov[i] = new NucleotideChannelCovariances();
        }
        
        for (int index = 0; index < recoveredICs.size(); index++) {
            NucleotideSequence control = controls.get(icids.get(index).intValue());
            FourIntensity[] intensities = recoveredICs.get(index);
            int seqlength = (control.size() < intensities.length) ? control.size() : intensities.length;

            for (int cycle = 0; cycle < seqlength; cycle++) {
                Nucleotide nuc = control.get(cycle);
                FourIntensity sig = intensities[cycle];
                NucleotideChannelMeans mus = cmeans[cycle];

                ccov[cycle].add(nuc, sig, mus);
            }
        }

        System.out.println("Found " + recoveredICs.size() + " in " + numReads + " reads.");

        // Extend IC solutions to end of read length (if ICs are shorter than read)
        System.err.println("Extending cycle " + (controls.get(0).size() - 1) + " parameters from cycle " + controls.get(0).size() + " to cycle " + recoveredICs.get(0).length);
        for (int cycle = controls.get(0).size(); cycle < recoveredICs.get(0).length; cycle++) {
            cmeans[cycle] = cmeans[controls.get(0).size() - 1];
            ccov[cycle] = ccov[controls.get(0).size() - 1];
        }

        // Now compute probabilities for the bases
        Algebra alg = new Algebra();

        for (int cycle = 0; cycle < recoveredICs.get(0).length; cycle++) {
            ccov[cycle].invert();
        }

        FirecrestFileParser ffp2 = new FirecrestFileParser(FIRECREST_DIR, LANE, CYCLE_START, CYCLE_STOP);

        SAMFileHeader sfh = new SAMFileHeader();
        SAMFileWriter sfw = new SAMFileWriterFactory().makeSAMOrBAMWriter(sfh, false, SAM_OUT);

        numReads = 0;
        while (ffp2.hasNext()) {
            if (numReads % 1000000 == 0) {
                System.err.println("Basecalled " + numReads + " reads.");
            }

            FourIntensity[] intensities = ffp2.next().getIntensities();

            byte[] asciiseq = new byte[intensities.length];
            byte[] bestqual = new byte[intensities.length];
            byte[] nextbestqual = new byte[intensities.length];

            for (int cycle = 0; cycle < intensities.length; cycle++) {
                FourIntensity fi = intensities[cycle];

                double[] likes = new double[4];
                double total = 0.0;

                for (Nucleotide nuc : Nucleotide.values()) {
                    double norm = Math.sqrt(alg.det(ccov[cycle].channelCovariances(nuc)))/Math.pow(2.0*Math.PI, 2.0);

                    DoubleMatrix1D sub = subtract(fi, cmeans[cycle].channelMeans(nuc));
                    DoubleMatrix1D Ax = alg.mult(ccov[cycle].channelCovariances(nuc), sub);

                    double exparg = -0.5*alg.mult(sub, Ax);
                    likes[nuc.ordinal()] = norm*Math.exp(exparg);
                    total += likes[nuc.ordinal()];
                }

                Nucleotide call1 = Nucleotide.A;
                double prob1 = likes[0]/total;
                for (int i = 1; i < 4; i++) {
                    if (likes[i]/total > prob1) {
                        prob1 = likes[i]/total;

                        switch (i) {
                            case 1: call1 = Nucleotide.C; break;
                            case 2: call1 = Nucleotide.G; break;
                            case 3: call1 = Nucleotide.T; break;
                        }
                    }
                }

                Nucleotide call2 = Nucleotide.A;
                double prob2 = 0.0;
                for (int i = 0; i < 4; i++) {
                    if (i != call1.ordinal() && likes[i]/total > prob2 && likes[i]/total < prob1) {
                        prob2 = likes[i]/total;

                        switch (i) {
                            case 0: call2 = Nucleotide.A; break;
                            case 1: call2 = Nucleotide.C; break;
                            case 2: call2 = Nucleotide.G; break;
                            case 3: call2 = Nucleotide.T; break;
                        }
                    }
                }

                asciiseq[cycle] = (byte) call1.asChar();
                bestqual[cycle] = toPhredScore(prob1);
                nextbestqual[cycle] = toCompressedQuality(call2, prob2);
            }

            SAMRecord sr = new SAMRecord(sfh);
            sr.setReadName(Integer.toString(numReads));
            sr.setReadUmappedFlag(true);
            sr.setReadBases(asciiseq);
            sr.setBaseQualities(bestqual);
            sr.setAttribute("SQ", nextbestqual);
            sfw.addAlignment(sr);

            numReads++;
        }

        sfw.close();

        System.out.println("Done.");
    }

    private static byte toPhredScore(double prob) {
        byte qual = (1.0 - prob < 0.00001) ? 40 : (byte) (-10*Math.log10(1.0 - prob));
        return (qual > 40) ? 40 : qual;
    }

    private static DoubleMatrix1D subtract(FourIntensity a, FourIntensity b) {
        DoubleMatrix1D sub = (DoubleFactory1D.dense).make(4);

        for (int i = 0; i < 4; i++) {
            sub.set(i, a.getChannelIntensity(i) - b.getChannelIntensity(i));
        }

        return sub;
    }

    private static byte toCompressedQuality(Nucleotide base, double prob) {
        byte compressedQual = (byte) base.ordinal();
        byte cprob = (byte) (100.0*prob);
        byte qualmask = (byte) 252;
        compressedQual += ((cprob << 2) & qualmask);

        return compressedQual;
    }

    private static NucleotideSequence toNucleotideSequence(FourIntensity[] intensities) {
        NucleotideSequence ns = new NucleotideSequence(intensities.length);

        for (int cycle = 0; cycle < intensities.length; cycle++) {
            int brightestChannel = intensities[cycle].brightestChannel();

            Nucleotide nt = Nucleotide.A;
            switch (brightestChannel) {
                case 0: nt = Nucleotide.A; break;
                case 1: nt = Nucleotide.C; break;
                case 2: nt = Nucleotide.G; break;
                case 3: nt = Nucleotide.T; break;
            }

            ns.set(cycle, nt);
        }

        return ns;
    }
}
