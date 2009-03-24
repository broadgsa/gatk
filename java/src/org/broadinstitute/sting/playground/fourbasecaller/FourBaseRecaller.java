package org.broadinstitute.sting.playground.fourbasecaller;

import java.io.File;
import java.io.FilenameFilter;
import java.io.FileFilter;
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
import edu.mit.broad.picard.illumina.BustardFileParser;
import edu.mit.broad.picard.illumina.BustardReadData;
import edu.mit.broad.picard.illumina.BustardFileParser_1_1;

public class FourBaseRecaller {
    public static void main(String[] argv) {
        // Parse args
        File FIRECREST_DIR = new File(argv[0]);
        int LANE = Integer.valueOf(argv[1]);
        File SAM_OUT = new File(argv[2]);
        int CYCLE_START = Integer.valueOf(argv[3]);
        int CYCLE_STOP = Integer.valueOf(argv[4]);
        boolean isPaired = Boolean.valueOf(argv[5]);

        int readLength = (CYCLE_STOP - CYCLE_START);
        File BUSTARD_DIR = getBustardDirectory(FIRECREST_DIR);
        int limit = 1000000;

        NucleotideChannelMeans[] cmeans = new NucleotideChannelMeans[readLength];
        NucleotideChannelCovariances[] ccov = new NucleotideChannelCovariances[readLength];
        for (int i = 0; i < readLength; i++) {
            cmeans[i] = new NucleotideChannelMeans();
            ccov[i] = new NucleotideChannelCovariances();
        }

        // Loop through bustard data and compute signal means
        FirecrestFileParser ffp1 = new FirecrestFileParser(FIRECREST_DIR, LANE, CYCLE_START, CYCLE_STOP);
        BustardFileParser_1_1 bfp1 = new BustardFileParser_1_1(BUSTARD_DIR, LANE, isPaired, "BS");

        for (int queryid = 0; queryid < limit && ffp1.hasNext(); queryid++) {
            if (queryid % (limit/10) == 0) {
                System.err.println("Processed " + queryid + " reads for means.");
            }
            FourIntensity[] intensities = ffp1.next().getIntensities();
            String rsq = (CYCLE_START == 0) ? bfp1.next().getFirstReadSequence() : bfp1.next().getSecondReadSequence();

            for (int cycle = 0; cycle < readLength; cycle++) {
                FourIntensity sig = intensities[cycle];

                if      (rsq.charAt(cycle) == 'A') { cmeans[cycle].add(Nucleotide.A, sig); }
                else if (rsq.charAt(cycle) == 'C') { cmeans[cycle].add(Nucleotide.C, sig); }
                else if (rsq.charAt(cycle) == 'G') { cmeans[cycle].add(Nucleotide.G, sig); }
                else if (rsq.charAt(cycle) == 'T') { cmeans[cycle].add(Nucleotide.T, sig); }
            }
        }

        // Go through the data again and compute signal covariances
        FirecrestFileParser ffp2 = new FirecrestFileParser(FIRECREST_DIR, LANE, CYCLE_START, CYCLE_STOP);
        BustardFileParser_1_1 bfp2 = new BustardFileParser_1_1(BUSTARD_DIR, LANE, isPaired, "BS");

        for (int queryid = 0; queryid < limit && ffp2.hasNext(); queryid++) {
            if (queryid % (limit/10) == 0) {
                System.err.println("Processed " + queryid + " reads for covariances.");
            }
            
            FourIntensity[] intensities = ffp2.next().getIntensities();
            String rsq = (CYCLE_START == 0) ? bfp2.next().getFirstReadSequence() : bfp2.next().getSecondReadSequence();

            for (int cycle = 0; cycle < readLength; cycle++) {
                FourIntensity sig = intensities[cycle];
                NucleotideChannelMeans mus = cmeans[cycle];

                if      (rsq.charAt(cycle) == 'A') { ccov[cycle].add(Nucleotide.A, sig, mus); }
                else if (rsq.charAt(cycle) == 'C') { ccov[cycle].add(Nucleotide.C, sig, mus); }
                else if (rsq.charAt(cycle) == 'G') { ccov[cycle].add(Nucleotide.G, sig, mus); }
                else if (rsq.charAt(cycle) == 'T') { ccov[cycle].add(Nucleotide.T, sig, mus); }
            }
        }

        // Now compute probabilities for the bases
        Algebra alg = new Algebra();

        for (int cycle = 0; cycle < readLength; cycle++) {
            ccov[cycle].invert();
        }

        FirecrestFileParser ffp3 = new FirecrestFileParser(FIRECREST_DIR, LANE, CYCLE_START, CYCLE_STOP);

        SAMFileHeader sfh = new SAMFileHeader();
        SAMFileWriter sfw = new SAMFileWriterFactory().makeSAMOrBAMWriter(sfh, false, SAM_OUT);

        for (int queryid = 0; ffp3.hasNext(); queryid++) {
            if (queryid % limit == 0) {
                System.err.println("Basecalled " + queryid + " reads.");
            }

            FourIntensity[] intensities = ffp3.next().getIntensities();

            byte[] asciiseq = new byte[readLength];
            byte[] bestqual = new byte[readLength];
            byte[] nextbestqual = new byte[readLength];

            for (int cycle = 0; cycle < readLength; cycle++) {
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
            sr.setReadName(Integer.toString(queryid));
            sr.setReadUmappedFlag(true);
            sr.setReadBases(asciiseq);
            sr.setBaseQualities(bestqual);
            sr.setAttribute("SQ", nextbestqual);
            sfw.addAlignment(sr);

            queryid++;
        }

        sfw.close();

        System.err.println("Done.");
    }

    private static byte toPhredScore(double prob) {
        byte qual = (1.0 - prob < 0.00001) ? 40 : (byte) (-10*Math.log10(1.0 - prob));
        //System.out.println("prob=" + prob + " qual=" + qual);
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

    private static File getBustardDirectory(File firecrestDir) {
        FileFilter filter = new FileFilter() {
            public boolean accept(File file) {
                return (file.isDirectory() && file.getName().contains("Bustard"));
            }
        };

        File[] bustardDirs = firecrestDir.listFiles(filter);

        return bustardDirs[0];
    }
}
