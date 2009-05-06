package org.broadinstitute.sting.playground.fourbasecaller;

import edu.mit.broad.picard.illumina.BustardFileParser;
import edu.mit.broad.picard.illumina.BustardReadData;
import edu.mit.broad.picard.illumina.AbstractBustardFileParser;
import edu.mit.broad.picard.illumina.BustardFileParser_1_1;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.playground.illumina.FirecrestFileParser;
import org.broadinstitute.sting.playground.illumina.FirecrestReadData;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

public class FourBaseRecaller extends CommandLineProgram {
    public static FourBaseRecaller Instance = null;
    
    @Argument(fullName="dir",         shortName="D", doc="Illumina Bustard directory")
    public File DIR;
    @Argument(fullName="lane",        shortName="L", doc="Illumina flowcell lane")
    public int LANE;
    @Argument(fullName="run_barcode", shortName="B", doc="Illumina Run Barcode (e.g. 305PJAAXX080716)")
    public String RUN_BARCODE;
    @Argument(fullName="out",         shortName="O", doc="Output path for sam file")
    public File OUT;
    @Argument(fullName="end",         shortName="E", doc="End of read to process (0 = whole read, i.e. unpaired; 1 = first end; 2 = second end)", required=false)
    public int END = 0;
    @Argument(fullName="tlim",        shortName="T", doc="Number of reads to use for parameter initialization", required=false)
    public int TRAINING_LIMIT = 1000000000;
    @Argument(fullName="clim",        shortName="C", doc="Number of reads to basecall", required=false)
    public int CALLING_LIMIT = 1000000000;
    @Argument(fullName="raw",         shortName="R", doc="Use raw intensities?", required=false)
    public Boolean RAW = false;
    @Argument(fullName="old",         shortName="1", doc="Old Bustard 1.1 mode?", required=false)
    public Boolean OLD = false;
    @Argument(fullName="context",     shortName="X", doc="Correct for context?", required=false)
    public Boolean CONTEXT = false;

    public static void main(String[] argv) {
        Instance = new FourBaseRecaller();
        start(Instance, argv);
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

        AbstractBustardFileParser bfp;
        BustardReadData bread;

        FirecrestFileParser ffp;
        FirecrestReadData fread;

        bfp = OLD ? new BustardFileParser_1_1(DIR, LANE, isPaired, RUN_BARCODE) : new BustardFileParser(DIR, LANE, isPaired, RUN_BARCODE);
        bread = bfp.next();
        ffp = new FirecrestFileParser(DIR.getParentFile(), LANE);
        fread = ffp.next();

        int cycle_offset = (END <= 1) ? 0 : bread.getIntensities().length/2;
        BasecallingReadModel model = new BasecallingReadModel(bread.getFirstReadSequence().length(), CONTEXT);
        int queryid;

        // learn mean parameters
        System.err.println("Computing mean parameters...");

        queryid = 0;
        do {
            assert(bread.getXCoordinate() == fread.getXCoordinate() && bread.getYCoordinate() == fread.getYCoordinate());

            String bases = (END <= 1) ? bread.getFirstReadSequence() : bread.getSecondReadSequence();
            byte[] quals = (END <= 1) ? bread.getFirstReadPhredBinaryQualities() : bread.getSecondReadPhredBinaryQualities();
            double[][] intensities = bread.getIntensities();
            double[][] rawintensities = fread.getIntensities();

            for (int cycle = 0; cycle < bases.length(); cycle++) {
                char baseCur  = bases.charAt(cycle);
                byte qualCur  = quals[cycle];
                //double[] fourprob = getBaseProbabilityDistribution(bases.charAt(cycle), quals[cycle]);
                double[] fourintensity = (RAW || OLD) ? rawintensities[cycle + cycle_offset] : intensities[cycle + cycle_offset];

                model.addMeanPoint(cycle, cycle == 0 ? '.' : bases.charAt(cycle - 1), baseCur, qualCur, fourintensity);
            }

            queryid++;
        } while (queryid < TRAINING_LIMIT && bfp.hasNext() && (bread = bfp.next()) != null && (fread = ffp.next()) != null);

        // learn covariance parameters
        System.err.println("Computing covariance parameters...");

        bfp = OLD ? new BustardFileParser_1_1(DIR, LANE, isPaired, RUN_BARCODE) : new BustardFileParser(DIR, LANE, isPaired, RUN_BARCODE);
        bread = bfp.next();
        ffp = new FirecrestFileParser(DIR.getParentFile(), LANE);
        fread = ffp.next();

        queryid = 0;
        do {
            assert(bread.getXCoordinate() == fread.getXCoordinate() && bread.getYCoordinate() == fread.getYCoordinate());

            String bases = (END <= 1) ? bread.getFirstReadSequence() : bread.getSecondReadSequence();
            byte[] quals = (END <= 1) ? bread.getFirstReadPhredBinaryQualities() : bread.getSecondReadPhredBinaryQualities();
            double[][] intensities = bread.getIntensities();
            double[][] rawintensities = fread.getIntensities();

            for (int cycle = 0; cycle < bases.length(); cycle++) {
                char baseCur  = bases.charAt(cycle);
                byte qualCur  = quals[cycle];
                //double[] fourprob = getBaseProbabilityDistribution(bases.charAt(cycle), quals[cycle]);
                double[] fourintensity = (RAW || OLD) ? rawintensities[cycle + cycle_offset] : intensities[cycle + cycle_offset];

                model.addCovariancePoint(cycle, cycle == 0 ? '.' : bases.charAt(cycle - 1), baseCur, qualCur, fourintensity);
            }

            queryid++;
        } while (queryid < TRAINING_LIMIT && bfp.hasNext() && (bread = bfp.next()) != null && (fread = ffp.next()) != null);

        // write model parameters
        model.write(debugdir);

        // call bases
        System.err.println("Calling bases...");

        SAMFileHeader sfh = new SAMFileHeader();
        SAMFileWriter sfw = new SAMFileWriterFactory().makeSAMOrBAMWriter(sfh, false, OUT);
        
        bfp = OLD ? new BustardFileParser_1_1(DIR, LANE, isPaired, RUN_BARCODE) : new BustardFileParser(DIR, LANE, isPaired, RUN_BARCODE);
        bread = bfp.next();
        ffp = new FirecrestFileParser(DIR.getParentFile(), LANE);
        fread = ffp.next();

        int[][] base_counts = new int[4][bread.getFirstReadSequence().length()];
        int bases_incorrect = 0, bases_total = 0;

        if (debugout != null) { debugout.println("cycle int_a int_c int_g int_t bustard_base kiran_base bustard_prob kiran_prob kiran_base_prev"); }

        queryid = 0;
        do {
            assert(bread.getXCoordinate() == fread.getXCoordinate() && bread.getYCoordinate() == fread.getYCoordinate());

            String bases = (END <= 1) ? bread.getFirstReadSequence() : bread.getSecondReadSequence();
            byte[] quals = (END <= 1) ? bread.getFirstReadPhredBinaryQualities() : bread.getSecondReadPhredBinaryQualities();
            double[][] intensities = bread.getIntensities();
            double[][] rawintensities = fread.getIntensities();

            byte[] asciiseq = new byte[bases.length()];
            byte[] bestqual = new byte[bases.length()];
            byte[] nextbestqual = new byte[bases.length()];

            for (int cycle = 0; cycle < bases.length(); cycle++) {
                double[] fourintensity = (RAW || OLD) ? rawintensities[cycle + cycle_offset] : intensities[cycle + cycle_offset];

                char basePrev = (cycle == 0 ? '.' : (char) asciiseq[cycle - 1]);
                byte qualPrev = (cycle == 0 ? 40 : bestqual[cycle - 1]);

                FourProb fp = model.computeProbabilities(cycle, basePrev, qualPrev, fourintensity);

                asciiseq[cycle] = (byte) fp.baseAtRank(0);
                bestqual[cycle] = fp.qualAtRank(0);
                nextbestqual[cycle] = QualityUtils.baseAndProbToCompressedQuality(fp.indexAtRank(1), fp.probAtRank(1));

                if (debugout != null && bases.charAt(cycle) != '.' && base_counts[BaseUtils.simpleBaseToBaseIndex(bases.charAt(cycle))][cycle] < 1000) {
                    debugout.println(cycle + " " + fourintensity[0] + " " + fourintensity[1] + " " + fourintensity[2] + " " + fourintensity[3] + " " + (bases.charAt(cycle)) + " " + ((char) asciiseq[cycle]) + " " + bestqual[cycle] + " " + quals[cycle] + " " + basePrev);

                    base_counts[BaseUtils.simpleBaseToBaseIndex(bases.charAt(cycle))][cycle]++;
                }

                bases_incorrect += (bases.charAt(cycle) == (char) asciiseq[cycle]) ? 0 : 1;
                bases_total++;
            }

            sfw.addAlignment(constructSAMRecord("KIR_", new String(asciiseq), bestqual, nextbestqual, bread, sfh));
            sfw.addAlignment(constructSAMRecord("BUS_", bases,                quals,    null,         bread, sfh));

            queryid++;
        } while (queryid < CALLING_LIMIT && bfp.hasNext() && (bread = bfp.next()) != null && (fread = ffp.next()) != null);

        if (debugout != null) {
            debugout.close();
        }
        sfw.close();

        System.out.println("Disagreement rate: " + ((double) bases_incorrect)/((double) bases_total));

        return 0;
    }

    private double[] getBaseProbabilityDistribution(char base, byte qual) {
        double[] dist = new double[4];

        int baseIndex = BaseUtils.simpleBaseToBaseIndex(base);
        dist[baseIndex] = QualityUtils.qualToProb(qual);
        double residualprob = (1.0 - dist[baseIndex])/3.0;

        for (int i = 0; i < 4; i++) {
            if (i != baseIndex) {
                dist[i] = residualprob;
            }
        }

        return dist;
    }

    private SAMRecord constructSAMRecord(String readNamePrefix, String bases, byte[] bestqual, byte[] nextbestqual, BustardReadData bread, SAMFileHeader sfh) {
        SAMRecord sr = new SAMRecord(sfh);
        
        sr.setReadName(readNamePrefix + bread.getReadName());
        sr.setReadUmappedFlag(true);
        sr.setReadString(bases);
        sr.setBaseQualities(bestqual);
        if (nextbestqual != null) { sr.setAttribute("SQ", nextbestqual); }
        sr.setReadFailsVendorQualityCheckFlag(!bread.isPf());

        return sr;
    }
}
