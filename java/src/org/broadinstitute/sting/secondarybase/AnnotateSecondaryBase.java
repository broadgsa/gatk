package org.broadinstitute.sting.secondarybase;

import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.secondarybase.BasecallingReadModel;

import java.io.File;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Vector;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

public class AnnotateSecondaryBase extends CommandLineProgram {
    public static AnnotateSecondaryBase Instance = null;

    @Argument(fullName="dir", shortName="D", doc="Illumina Firecrest directory") public File FIRECREST_DIR;
    @Argument(fullName="lane", shortName="L", doc="Illumina flowcell lane") public int LANE;
    @Argument(fullName="sam_in", shortName="SI", doc="The file to annotate") public File SAM_IN;
    @Argument(fullName="sam_out", shortName="SO", doc="Output path for sam file") public File SAM_OUT;
    @Argument(fullName="reference", shortName="R", doc="Reference sequence to which sam_in is aligned (in fasta format)") public File REFERENCE;
    @Argument(fullName="cycle_begin", shortName="CB", doc="On what cycle does the read begin? (0-based inclusive)") public int CYCLE_BEGIN;
    @Argument(fullName="cycle_end", shortName="CE", doc="On what cycle does the read end? (0-based inclusive)") public int CYCLE_END;
    @Argument(fullName="tlim", shortName="T", doc="Number of reads to use for parameter initialization", required=false) public int TRAINING_LIMIT = 1000000;
    @Argument(fullName="clim", shortName="C", doc="Number of reads to basecall", required=false) public int CALLING_LIMIT = Integer.MAX_VALUE;
    @Argument(fullName="context", shortName="X", doc="Attempt to correct for context?", required=false) public Boolean CONTEXT = false;

    public static void main(String[] argv) {
        Instance = new AnnotateSecondaryBase();
        start(Instance, argv);
    }

    protected int execute() {
        // Iterate through bam file
            // take an alignment
            // verify that it has zero mismatches to the reference
            // if so, store it
            // continue until we've picked up TRAINING_LIMIT alignments
        Vector< HashMap<String, SAMRecord> > trainingSet = loadPreAlignedTrainingSet();

        // Iterate through the stored alignments
            // find the corresponding firecrest intensity line
            // add the info to the BasecallingReadModel
        BasecallingReadModel model = train(trainingSet);

        // Iterate through the bam file
            // take an alignment
            // call it
            // stick the 2nd-best bases into the file
            // spit the alignment back out

        return 0;        
    }

    private Vector< HashMap<String, SAMRecord> > loadPreAlignedTrainingSet() {
        FastaSequenceFile2 ref = new FastaSequenceFile2(REFERENCE);
        String currentContig = "none";
        byte[] refbases = null;

        SAMFileReader sf = new SAMFileReader(SAM_IN);
        SAMRecord sr;
        CloseableIterator<SAMRecord> sfit = sf.iterator();

        int numTrainingReads = 0;
        Vector< HashMap<String, SAMRecord> > trainingReads = new Vector< HashMap<String, SAMRecord> >(100);

        while (numTrainingReads < TRAINING_LIMIT && (sr = sfit.next()) != null) {
            if (sr.getCigar().numCigarElements() == 1) {
                int offset = sr.getAlignmentStart();

                if (!currentContig.matches(sr.getReferenceName())) {
                    ref.seekToContig(sr.getReferenceName());
                    currentContig = sr.getReferenceName();
                    refbases = ref.nextSequence().getBases();
                }

                int mismatches = 0;

                String refString = "";
                for (int i = offset, j = 0; i < offset + sr.getReadBases().length; i++, j++) {
                    refString += (char) refbases[i - 1];

                    mismatches += (BaseUtils.simpleBaseToBaseIndex((char) refbases[i - 1]) !=
                                   BaseUtils.simpleBaseToBaseIndex((char) sr.getReadBases()[j]))
                                  ? 1 : 0;
                }

                if (mismatches == 0) {
                    Pattern p = Pattern.compile(":(\\d):(\\d+):(\\d+):(\\d+)#");
                    Matcher m = p.matcher(sr.getReadName());

                    if (m.find()) {
                        int tile = Integer.valueOf(m.group(2));
                        String readKey = String.format("%s:%s:%s:%s", m.group(1), m.group(2), m.group(3), m.group(4));

                        if (tile > trainingReads.size()) {
                            trainingReads.setSize(tile + 1);
                        }

                        if (trainingReads.get(tile) == null) {
                            HashMap<String, SAMRecord> tileTrainingReads = new HashMap<String, SAMRecord>();
                            trainingReads.set(tile, tileTrainingReads);
                        }

                        trainingReads.get(tile).put(readKey, sr);
                        numTrainingReads++;
                    } else {
                        throw new StingException("Pattern '" + p.pattern() + "' does not match against read name '" + sr.getReadName() + "'");
                    }
                }
            }
        }

        System.out.println("Collected " + numTrainingReads + " reads.");

        return trainingReads;
    }

    private BasecallingReadModel train(Vector< HashMap<String, SAMRecord> > trainingSet) {
        BasecallingReadModel model = new BasecallingReadModel(CYCLE_END - CYCLE_BEGIN + 1, CONTEXT);

        return model;
    }
}
