package org.broadinstitute.sting.playground.piecemealannotator;

import net.sf.samtools.*;
import org.broadinstitute.sting.secondarybase.BasecallingReadModel;
import org.broadinstitute.sting.secondarybase.BasecallingStats;
import org.broadinstitute.sting.secondarybase.FourProbRead;
import org.broadinstitute.sting.secondarybase.RawRead;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;

import java.io.File;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class TileAnnotator extends CommandLineProgram {
    public static TileAnnotator instance = null;

    private String currentContig = "";
    private byte[] refbases;
    //private int lane;
    //private int tile;

    @Argument(fullName="sam_tile_in", shortName="STI", doc="SAM tile file", required=false) public File SAM_TILE_IN;
    @Argument(fullName="sam_tile_out", shortName="STO", doc="Annotated SAM tile output file") public File SAM_TILE_OUT;
    @Argument(fullName="reference", shortName="R", doc="The fasta reference") public File REFERENCE;
    @Argument(fullName="bustard_dir", shortName="D", doc="The Bustard directory") public File BUSTARD_DIR;
    @Argument(fullName="training_limit", shortName="TL", doc="Number of reads to train from", required=false) public int TRAINING_LIMIT = 10000;
    @Argument(fullName="run_barcode", shortName="RB", doc="Illumina run barcode") public String RUN_BARCODE;
    @Argument(fullName="cycle_ranges", shortName="CR", doc="Cycle ranges for single-end or paired reads (i.e. '0-75,76-151') (0-based, inclusive)") public String CYCLE_RANGES;
    @Argument(fullName="lane", shortName="L", doc="The lane to process (if not specified, this will be read from the 'sam_tile_in' file)", required=false) public Integer lane;
    @Argument(fullName="tile", shortName="T", doc="The tile to process (if not specified, this will be read from the 'sam_tile_in' file)", required=false) public Integer tile;

    public static void main(String[] argv) {
        instance = new TileAnnotator();
        start(instance, argv);
    }

    protected int execute() {
        ArrayList<Pair<Integer, Integer>> cycleRanges = getCycleRanges(CYCLE_RANGES);

        System.out.printf("%s: Loading training set...\n", (new Date()).toString());
        ArrayList<RawRead> trainingData = loadTrainingData();

        System.out.printf("%s: Training model...\n", (new Date()).toString());
        BasecallingReadModel model = new BasecallingReadModel(trainingData);

        System.out.printf("%s: Calling bases...\n", (new Date()).toString());
        //callBases(model, cycleRanges, trainingData);
        callBases(model, cycleRanges);

        System.out.println("Done.");

        return 0;
    }

    //private void callBases(BasecallingReadModel model, ArrayList<Pair<Integer, Integer>> cycleRanges, ArrayList<RawRead> trainingData) {
    private void callBases(BasecallingReadModel model, ArrayList<Pair<Integer, Integer>> cycleRanges) {
        SAMFileHeader sheader = new SAMFileHeader();
        sheader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        SAMFileWriter swriter =  new SAMFileWriterFactory().makeSAMOrBAMWriter(sheader, true, SAM_TILE_OUT);

        IlluminaTile tileParser = new IlluminaTile(BUSTARD_DIR, lane, tile);

        BasecallingStats bstats = new BasecallingStats();

        //for (RawRead rr : trainingData) {
        for (RawRead rr : tileParser) {
            FourProbRead fpr = model.call(rr);

            //System.out.println(rr.getSequenceAsString());
            //System.out.println(fpr.getPrimaryBaseSequence());

            for (int rangeIndex = 0; rangeIndex < cycleRanges.size(); rangeIndex++) {
                FourProbRead fprEnd = fpr.getSubset(cycleRanges.get(rangeIndex).getFirst(), cycleRanges.get(rangeIndex).getSecond());
                RawRead rrEnd = rr.getSubset(cycleRanges.get(rangeIndex).getFirst(), cycleRanges.get(rangeIndex).getSecond());

                //System.out.println(rangeIndex + " : " + cycleRanges.get(rangeIndex).getFirst() + " " + cycleRanges.get(rangeIndex).getSecond());
                //System.out.println(rrEnd.getSequenceAsString());
                //System.out.println(fprEnd.getPrimaryBaseSequence());

                SAMRecord sr = constructSAMRecord(rrEnd, fprEnd, sheader, cycleRanges.size() > 1, rangeIndex == 1);

                swriter.addAlignment(sr);
            }

            //System.out.println();

            bstats.update(rr, fpr);
            bstats.notifyOnInterval(10000);
        }

        bstats.notifyNow();

        tileParser.close();
        swriter.close();
    }

    private ArrayList<Pair<Integer, Integer>> getCycleRanges(String cycleRangesString) {
        ArrayList< Pair<Integer, Integer> > cycleRanges = new ArrayList< Pair<Integer, Integer> >();

        String[] pieces = cycleRangesString.split(",");

        Pattern p = Pattern.compile("(\\d+)-(\\d+)");

        for (String piece : pieces) {
            Matcher m = p.matcher(piece);

            if (m.find()) {
                Integer cycleStart = new Integer(m.group(1));
                Integer cycleStop = new Integer(m.group(2));

                cycleRanges.add(new Pair<Integer, Integer>(cycleStart, cycleStop));
            }
        }

        if (cycleRanges.size() == 0) {
            throw new StingException("At least one cycle range must be specified.");
        }

        if (cycleRanges.size() > 2) {
            throw new StingException(cycleRanges.size() + " specified, but we're unable to handle more than 2.");
        }

        return cycleRanges;
    }

    private SAMRecord constructSAMRecord(RawRead rr, FourProbRead fpr, SAMFileHeader sfh, boolean isPaired, boolean isSecondEndOfPair) {
        SAMRecord sr = new SAMRecord(sfh);

        sr.setReadName(String.format("%s:%d:%d:%d:%d#0", RUN_BARCODE, lane, tile, rr.getXCoordinate(), rr.getYCoordinate()));
        sr.setReadUmappedFlag(true);
        sr.setReadString(rr.getSequenceAsString());
        sr.setBaseQualities(rr.getQuals());
        sr.setAttribute("SQ", fpr.getSQTag(rr));

        sr.setReadPairedFlag(isPaired);
        if (isPaired) {
            sr.setMateUnmappedFlag(true);
            sr.setFirstOfPairFlag(!isSecondEndOfPair);
            sr.setSecondOfPairFlag(isSecondEndOfPair);
        }

        return sr;
    }

    private ArrayList<RawRead> loadTrainingData() {
        FastaSequenceFile2 ref = new FastaSequenceFile2(REFERENCE);
        HashMap<String, SAMRecord> srs = loadTileAlignments(ref);
        return loadGoodReads(srs, BUSTARD_DIR);
    }

    private HashMap<String, SAMRecord> loadTileAlignments(FastaSequenceFile2 ref) {
        HashMap<String, SAMRecord> srs = new HashMap<String, SAMRecord>();
        HashSet<String> seenEnds = new HashSet<String>();
        
        int numPerfect = 0;

        if (SAM_TILE_IN != null && SAM_TILE_IN.exists()) {
            SAMFileReader sreader = new SAMFileReader(SAM_TILE_IN);

            for (SAMRecord sr : sreader) {
                Pattern p = Pattern.compile(":(\\d+):(\\d+):(\\d+):(\\d+)#");
                Matcher m = p.matcher(sr.getReadName());

                if (m.find()) {
                    this.lane = Integer.valueOf(m.group(1));
                    this.tile = Integer.valueOf(m.group(2));
                    int x = Integer.valueOf(m.group(3));
                    int y = Integer.valueOf(m.group(4));
                    boolean end = sr.getReadPairedFlag() && sr.getSecondOfPairFlag();

                    String otherKey = String.format("%d:%d:%b", x, y, !end);
                    String currentKey = String.format("%d:%d:%b", x, y, end);

                    seenEnds.add(currentKey);

                    if (isWellAligned(sr, ref)) {
                        if (srs.containsKey(otherKey) || !seenEnds.contains(otherKey)) {
                            srs.put(currentKey, sr);
                        }

                        if (srs.containsKey(currentKey) && srs.containsKey(otherKey)) {
                            numPerfect++;
                            if (numPerfect % (TRAINING_LIMIT < 1000 ? TRAINING_LIMIT : 1000) == 0) {
                                System.out.println("  " + numPerfect + " well-aligned reads");
                            }
                        }
                    } else {
                        if (srs.containsKey(otherKey)) {
                            srs.remove(otherKey);
                        }
                    }
                }

                if (numPerfect >= TRAINING_LIMIT) { break; }
            }

            sreader.close();
        }

        return srs;
    }

    private boolean isWellAligned(SAMRecord sr, FastaSequenceFile2 ref) {
        boolean valid = false;
        int mismatches = 0;

        if (!sr.getReadUnmappedFlag() && sr.getCigar().numCigarElements() == 1) {
            if (!currentContig.matches(sr.getReferenceName())) {
                ref.seekToContig(sr.getReferenceName());
                currentContig = sr.getReferenceName();

                refbases = ref.nextSequence().getBases();
            }

            byte[] readbases = sr.getReadBases();
            int offset = sr.getAlignmentStart();

            if (offset + readbases.length < refbases.length) {
                valid = true;
    
                //String refString = "";
                for (int i = offset, j = 0; i < offset + readbases.length; i++, j++) {
                    //refString += (char) refbases[i - 1];

                    int refbase = BaseUtils.simpleBaseToBaseIndex((char) refbases[i - 1]);
                    int readbase = BaseUtils.simpleBaseToBaseIndex((char) readbases[j]);

                    mismatches += (refbase >= 0 && readbase >= 0 && refbase != readbase) ? 1 : 0;
                }
            }
        }

        return (valid && mismatches == 0);
    }

    private ArrayList<RawRead> loadGoodReads(HashMap<String, SAMRecord> srs, File bustardDir) {
        ArrayList<RawRead> trainingData = new ArrayList<RawRead>();
        BoundedScoringSet<RawRead> additionalData = new BoundedScoringSet<RawRead>(TRAINING_LIMIT);

        IlluminaTile tileParser = new IlluminaTile(bustardDir, lane, tile);

        int correlatedReads = 0;
        for (RawRead rr : tileParser) {
            String key1 = String.format("%d:%d:%b", rr.getXCoordinate(), rr.getYCoordinate(), false);
            String key2 = String.format("%d:%d:%b", rr.getXCoordinate(), rr.getYCoordinate(), true);

            if (srs.containsKey(key1) && srs.containsKey(key2)) {
                byte[] quals = new byte[rr.getReadLength()];
                for (int cycle = 0; cycle < rr.getReadLength(); cycle++) {
                    quals[cycle] = (byte) (BaseUtils.simpleBaseToBaseIndex((char) rr.getSequence()[cycle]) >= 0 ? 50 : 0);
                }
                rr.setQuals(quals);

                trainingData.add(rr);

                correlatedReads++;
                if (correlatedReads % (TRAINING_LIMIT < 1000 ? TRAINING_LIMIT : 1000) == 0) {
                    System.out.println("  " + correlatedReads + " intensity-correlated reads");
                }
            } else {
                additionalData.add(rr);
            }

            //trainingData.add(rr);

            //if (trainingData.size() % 10000 == 0) {
            //    System.out.printf("  loaded %d reads\n", trainingData.size());
            //}
        }
        
        //System.out.printf("  loaded %d reads\n", trainingData.size());

        //for (RawRead rr : additionalData.toArray(new RawRead[0])) {
        //    System.out.println(rr.getSequenceAsString());
        //}

        tileParser.close();

        System.out.printf("  found %d perfect reads with an optional reservoir of %d good reads\n", trainingData.size(), additionalData.size());

        RawRead[] qrs = additionalData.toArray(new RawRead[0]);
        int limit = (TRAINING_LIMIT - trainingData.size() < additionalData.size()) ? (TRAINING_LIMIT - trainingData.size()) : additionalData.size();
        for (int i = 0; i < limit; i++) {
            trainingData.add(qrs[i]);
        }

        return trainingData;
    }
}
