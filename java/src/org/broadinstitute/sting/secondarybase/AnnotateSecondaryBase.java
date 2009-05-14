package org.broadinstitute.sting.secondarybase;

import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.secondarybase.BasecallingReadModel;

import java.io.File;
import java.util.HashMap;
import java.util.Vector;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;

public class AnnotateSecondaryBase extends CommandLineProgram {
    public static AnnotateSecondaryBase Instance = null;

    @Argument(fullName="dir", shortName="D", doc="Illumina Bustard directory") public File BUSTARD_DIR;
    @Argument(fullName="lane", shortName="L", doc="Illumina flowcell lane") public int LANE;
    @Argument(fullName="sam_in", shortName="SI", doc="The file to use for training and annotation", required=false) public File SAM_IN;
    @Argument(fullName="sam_out", shortName="SO", doc="Output path for sam file") public File SAM_OUT;
    @Argument(fullName="reference", shortName="R", doc="Reference sequence to which sam_in is aligned (in fasta format)") public File REFERENCE;
    @Argument(fullName="cycle_begin", shortName="CB", doc="On what cycle does the read begin? (0-based inclusive)") public int CYCLE_BEGIN;
    @Argument(fullName="cycle_end", shortName="CE", doc="On what cycle does the read end? (0-based inclusive)") public int CYCLE_END;
    @Argument(fullName="tlim", shortName="T", doc="Number of reads to use for parameter initialization", required=false) public int TRAINING_LIMIT = 250000;
    @Argument(fullName="clim", shortName="C", doc="Number of reads to basecall", required=false) public int CALLING_LIMIT = Integer.MAX_VALUE;
    @Argument(fullName="context", shortName="X", doc="Attempt to correct for context?", required=false) public Boolean CONTEXT = false;
    @Argument(fullName="runbarcode", shortName="B", doc="Run barcode (embedded as part of the read name") public String RUN_BARCODE;

    public static void main(String[] argv) {
        Instance = new AnnotateSecondaryBase();
        start(Instance, argv);
    }

    protected int execute() {
        BasecallingTrainingSet trainingSet = new BasecallingTrainingSet(BUSTARD_DIR, LANE, CYCLE_BEGIN, CYCLE_END, TRAINING_LIMIT);

        if (SAM_IN == null || !SAM_IN.exists()) {
            // Iterate through raw Firecrest data and store the first N reads up to TRAINING_LIMIT
            System.out.println("Training from the first " + TRAINING_LIMIT + " reads in the raw data.");
            trainingSet.loadFirstNUnambiguousReadsTrainingSet();
        } else {
            // Find alignments with zero mismatches and store them until we've picked up TRAINING_LIMIT alignments
            System.out.println("Training from the first " + TRAINING_LIMIT + " perfect reads in the aligned data.");
            trainingSet.loadPreAlignedTrainingSet(SAM_IN, REFERENCE);
        }

        // Iterate through the stored training data and add the info to the BasecallingReadModel
        BasecallingReadModel model = new BasecallingReadModel(CYCLE_END - CYCLE_BEGIN + 1, CONTEXT);
        model.train(trainingSet);

        // Call bases and write results
        SAMFileHeader sfh = new SAMFileHeader();
        SAMFileWriter sfw = new SAMFileWriterFactory().makeSAMOrBAMWriter(sfh, false, SAM_OUT);

        IlluminaParser iparser = new IlluminaParser(BUSTARD_DIR, LANE, CYCLE_BEGIN, CYCLE_END);

        RawRead rr;
        while ((rr = iparser.next()) != null) {
            FourProbRead fpr = model.call(rr);

            System.out.println(rr.getSequenceAsString());
            System.out.println(fpr.getPrimaryBaseSequence());

            SAMRecord sr = constructSAMRecord(rr, fpr, sfh, RUN_BARCODE);
            sfw.addAlignment(sr);
        }

        iparser.close();
        
        sfw.close();

        return 0;        
    }

    private SAMRecord constructSAMRecord(RawRead rr, FourProbRead fpr, SAMFileHeader sfh, String runBarcode) {
        SAMRecord sr = new SAMRecord(sfh);

        sr.setReadName(runBarcode + ":" + rr.getReadKey() + "#0");
        sr.setReadUmappedFlag(true);
        sr.setReadString(rr.getSequenceAsString());
        sr.setBaseQualities(rr.getQuals());
        sr.setAttribute("SQ", fpr.getSQTag(rr));

        return sr;
    }
}
