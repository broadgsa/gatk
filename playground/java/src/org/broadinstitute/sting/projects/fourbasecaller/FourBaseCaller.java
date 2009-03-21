package org.broadinstitute.sting.projects.fourbasecaller;

import java.io.File;
import java.util.Vector;
import java.lang.Math;

import org.broadinstitute.sting.illumina.FirecrestFileParser;
import org.broadinstitute.sting.illumina.FourIntensity;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleFactory1D;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;

public class FourBaseCaller {
    public static void main(String[] args) {
        // Parse args
        File firecrestDirectory = new File(args[0]);
        int lane = Integer.valueOf(args[1]);
        String icrefpath = args[2];

        // Load ic reference
        ManyNucleotideSequences controls = new ManyNucleotideSequences(icrefpath);
        controls.reverseComplement();

        // Loop through firecrest data, find ICs, and compute signal means
        Vector<FourIntensity[]> recoveredICs = new Vector<FourIntensity[]>();
        Vector<Integer> icids = new Vector<Integer>();

        FirecrestFileParser ffp = new FirecrestFileParser(firecrestDirectory, lane);

        NucleotideChannelMeans[] cmeans = new NucleotideChannelMeans[controls.get(0).size()];
        for (int i = 0; i < cmeans.length; i++) {
            cmeans[i] = new NucleotideChannelMeans();
        }

        int numReads = 0;
        while (ffp.hasNext() && numReads < 100000) {
            if (numReads % 10000 == 0) {
                System.out.println("Processed " + numReads + " reads.");
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

        // Now compute probabilities for the bases
        Algebra alg = new Algebra();

        for (int cycle = 0; cycle < recoveredICs.get(0).length; cycle++) {
            ccov[cycle].invert();
        }

        FirecrestFileParser ffp2 = new FirecrestFileParser(firecrestDirectory, lane);

        numReads = 0;

        SAMFileHeader sfh = new SAMFileHeader();
        SAMFileWriter sfw = new SAMFileWriterFactory().makeSAMOrBAMWriter(sfh, false, new File("/wga/scr1/GSA/kiran/basecalling/simulation/test.sam"));

        while (ffp2.hasNext() && numReads < 100000) {
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

                //System.out.print(call.asChar());
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
    }

    private static byte toPhredScore(double prob) {
        byte qual = (prob - 1.0 < 0.00001) ? 40 : (byte) (-10*Math.log10(1.0 - prob));
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
