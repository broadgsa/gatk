package org.broadinstitute.sting.playground.fourbasecaller;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;
import org.broadinstitute.sting.utils.QualityUtils;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import edu.mit.broad.picard.illumina.BustardFileParser;
import edu.mit.broad.picard.illumina.BustardReadData;

public class FourBaseRecaller extends CommandLineProgram {
    public static FourBaseRecaller Instance = null;
    
    public File DIR;
    public int LANE;
    public File OUT;
    public int END = 0;
    public int TRAINING_LIMIT = 1000000000;
    public int CALLING_LIMIT = 1000000000;

    public static void main(String[] argv) {
        Instance = new FourBaseRecaller();
        start(Instance, argv);
    }

    protected void setupArgs() {
        m_parser.addRequiredArg("dir",  "D", "Illumina Bustard directory", "DIR");
        m_parser.addRequiredArg("lane", "L", "Illumina flowcell lane", "LANE");
        m_parser.addRequiredArg("out",  "O", "Output path for sam file", "OUT");
        m_parser.addOptionalArg("end",  "E", "End of read to process (0 = whole read, i.e. unpaired; 1 = first end; 2 = second end)", "END");
        m_parser.addOptionalArg("tlim", "T", "Number of reads to use for parameter initialization", "TRAINING_LIMIT");
        m_parser.addOptionalArg("clim", "C", "Number of reads to basecall", "CALLING_LIMIT");
    }

    protected int execute() {
        boolean isPaired = (END > 0);

        // Set up debugging paths
        File debugdir = new File(OUT.getPath() + ".debug/");
        debugdir.mkdir();
        PrintWriter debugout = null;
        try {
            debugout = new PrintWriter(debugdir.getPath() + "/debug.out");
        } catch (IOException e) {
        }

        BustardFileParser bfp;
        BustardReadData bread;

        bfp = new BustardFileParser(DIR, LANE, isPaired, "FB");
        bread = bfp.next();

        int cycle_offset = (END <= 1) ? 0 : bread.getIntensities().length/2;
        BasecallingReadModel model = new BasecallingReadModel(bread.getFirstReadSequence().length());
        int queryid;

        // learn mean parameters
        if (debugout != null) { debugout.println("intensity int_a int_c int_g int_t base"); }

        queryid = 0;
        do {
            String bases = (END <= 1) ? bread.getFirstReadSequence() : bread.getSecondReadSequence();
            byte[] quals = (END <= 1) ? bread.getFirstReadPhredBinaryQualities() : bread.getSecondReadPhredBinaryQualities();
            double[][] intensities = bread.getIntensities();

            for (int cycle = 0; cycle < bases.length(); cycle++) {
                char baseCur  = bases.charAt(cycle);
                byte qualCur  = quals[cycle];
                double[] fourintensity = intensities[cycle + cycle_offset];

                if (debugout != null && cycle == 0) {
                    debugout.println("intensity " + intensities[0][0] + " " + intensities[0][1] + " " + intensities[0][2] + " " + intensities[0][3] + " " + baseCur);
                }

                model.addMeanPoint(cycle, baseCur, qualCur, fourintensity);
            }

            queryid++;
        } while (queryid < TRAINING_LIMIT && bfp.hasNext() && (bread = bfp.next()) != null);

        // learn covariance parameters
        bfp = new BustardFileParser(DIR, LANE, isPaired, "FB");
        bread = bfp.next();

        queryid = 0;
        do {
            String bases = (END <= 1) ? bread.getFirstReadSequence() : bread.getSecondReadSequence();
            byte[] quals = (END <= 1) ? bread.getFirstReadPhredBinaryQualities() : bread.getSecondReadPhredBinaryQualities();
            double[][] intensities = bread.getIntensities();

            for (int cycle = 0; cycle < bases.length(); cycle++) {
                char baseCur  = bases.charAt(cycle);
                byte qualCur  = quals[cycle];
                double[] fourintensity = intensities[cycle + cycle_offset];

                model.addCovariancePoint(cycle, baseCur, qualCur, fourintensity);
            }

            queryid++;
        } while (queryid < TRAINING_LIMIT && bfp.hasNext() && (bread = bfp.next()) != null);

        // write debugging info
        model.write(debugdir);

        // call bases
        SAMFileHeader sfh = new SAMFileHeader();
        SAMFileWriter sfw = new SAMFileWriterFactory().makeSAMOrBAMWriter(sfh, false, OUT);
        
        bfp = new BustardFileParser(DIR, LANE, isPaired, "FB");
        bread = bfp.next();

        queryid = 0;
        do {
            String bases = (END <= 1) ? bread.getFirstReadSequence() : bread.getSecondReadSequence();
            byte[] quals = (END <= 1) ? bread.getFirstReadPhredBinaryQualities() : bread.getSecondReadPhredBinaryQualities();
            double[][] intensities = bread.getIntensities();

            byte[] asciiseq = new byte[bases.length()];
            byte[] bestqual = new byte[bases.length()];
            byte[] nextbestqual = new byte[bases.length()];

            for (int cycle = 0; cycle < bases.length(); cycle++) {
                double[] fourintensity = intensities[cycle + cycle_offset];

                FourProb fp = model.computeProbabilities(cycle, fourintensity);

                //if (cycle == 0) {
                //    System.out.println("result " + intensities[0][0] + " " + intensities[0][1] + " " + intensities[0][2] + " " + intensities[0][3] + " " + bases.charAt(0) + " " + fp.toString());
                //}

                asciiseq[cycle] = (byte) fp.baseAtRank(0);
                bestqual[cycle] = fp.qualAtRank(0);
                nextbestqual[cycle] = QualityUtils.baseAndProbToCompressedQuality(fp.indexAtRank(1), fp.probAtRank(1));
            }

            sfw.addAlignment(constructSAMRecord("KIR_", new String(asciiseq), bestqual, nextbestqual, isPaired, END, bread, sfh));
            sfw.addAlignment(constructSAMRecord("BUS_", bases,                quals,    null,         isPaired, END, bread, sfh));

            queryid++;
        } while (queryid < CALLING_LIMIT && bfp.hasNext() && (bread = bfp.next()) != null);

        sfw.close();

        return 0;
    }

    private SAMRecord constructSAMRecord(String readNamePrefix, String bases, byte[] bestqual, byte[] nextbestqual, boolean isPaired, int END, BustardReadData bread, SAMFileHeader sfh) {
        SAMRecord sr = new SAMRecord(sfh);
        
        sr.setReadName(readNamePrefix + bread.getReadName());
        sr.setReadUmappedFlag(true);
        sr.setReadString(bases);
        sr.setBaseQualities(bestqual);
        if (nextbestqual != null) { sr.setAttribute("SQ", nextbestqual); }
        sr.setReadFailsVendorQualityCheckFlag(!bread.isPf());
        if (isPaired) {
            sr.setMateUnmappedFlag(true);
            sr.setFirstOfPairFlag(END <= 1);
            sr.setFirstOfPairFlag(END > 1);
        }

        return sr;
    }
}
