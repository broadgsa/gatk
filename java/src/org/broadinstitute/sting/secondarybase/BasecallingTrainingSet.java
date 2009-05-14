package org.broadinstitute.sting.secondarybase;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.util.CloseableIterator;

import java.util.HashMap;
import java.util.Vector;
import java.util.ArrayList;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.io.File;

import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.StingException;

public class BasecallingTrainingSet {
    private File bustardDir;
    private int lane;
    private int cycleBegin;
    private int cycleEnd;
    private int trainingLimit;

    private ArrayList<RawRead> trainingData;

    public BasecallingTrainingSet(File bustardDir, int lane, int cycleBegin, int cycleEnd, int trainingLimit) {
        this.bustardDir = bustardDir;
        this.lane = lane;
        this.cycleBegin = cycleBegin;
        this.cycleEnd = cycleEnd;
        this.trainingLimit = trainingLimit;
    }

    public ArrayList<RawRead> getTrainingData() {
        return this.trainingData;
    }

    public void setTrainingData(ArrayList<RawRead> trainingData) {
        this.trainingData = trainingData;
    }

    public void loadFirstNUnambiguousReadsTrainingSet() {
        this.trainingData = new ArrayList<RawRead>(trainingLimit);

        IlluminaParser iparser = new IlluminaParser(bustardDir, lane, cycleBegin, cycleEnd);

        RawRead rawread;
        int numreads = 0;

        while (numreads < trainingLimit && (rawread = iparser.next()) != null) {
            int numAmbiguous = 0;
            byte[] sequence = rawread.getSequence();
            
            for ( byte byteBase : sequence ) {
                if (BaseUtils.simpleBaseToBaseIndex((char) byteBase) == -1) {
                    numAmbiguous++;
                }
            }

            if (numAmbiguous == 0) {
                trainingData.add(rawread);
                numreads++;
            }
        }
    }

    public void loadPreAlignedTrainingSet(File samIn, File reference) {
        Vector< HashMap<String, SAMRecord> > trainingReads = getPerfectAlignmentsByTile(samIn, reference);

        trainingData = correlateReadsAndIntensities(trainingReads);
    }

    private Vector<HashMap<String, SAMRecord>> getPerfectAlignmentsByTile(File samIn, File reference) {
        FastaSequenceFile2 ref = new FastaSequenceFile2(reference);
        String currentContig = "none";
        byte[] refbases = null;

        SAMFileReader sf = new SAMFileReader(samIn);
        SAMRecord sr;
        CloseableIterator<SAMRecord> sfit = sf.iterator();

        Vector< HashMap<String, SAMRecord> > trainingReads = new Vector< HashMap<String, SAMRecord> >(101);
        int numTrainingReads = 0;

        while (numTrainingReads < trainingLimit && (sr = sfit.next()) != null) {
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
                            trainingReads.set(tile, new HashMap<String, SAMRecord>());
                        }

                        trainingReads.get(tile).put(readKey, sr);
                        numTrainingReads++;
                    } else {
                        throw new StingException("Pattern '" + p.pattern() + "' does not match against read name '" + sr.getReadName() + "'");
                    }
                }
            }
        }

        return trainingReads;
    }

    private ArrayList<RawRead> correlateReadsAndIntensities(Vector<HashMap<String, SAMRecord>> trainingReads) {
        ArrayList<RawRead> newTrainingData = new ArrayList<RawRead>(trainingLimit);

        IlluminaParser iparser = new IlluminaParser(bustardDir, lane, cycleBegin, cycleEnd);

        for (int tile = 1; tile < trainingReads.size(); tile++) {
            iparser.seekToTile(tile);

            int tileReadCount = 0;

            RawRead iread;
            while (trainingReads.get(tile) != null && tileReadCount < trainingReads.get(tile).size() && (iread = iparser.next()) != null) {
                String readKey = iread.getReadKey();

                if (trainingReads.get(tile).containsKey(readKey)) {
                    byte[] quals = new byte[iread.getReadLength()];
                    for (int qualIndex = 0; qualIndex < quals.length; qualIndex++) {
                        quals[qualIndex] = 40;
                    }

                    iread.setQuals(quals);
                    newTrainingData.add(iread);

                    tileReadCount++;
                }
            }
        }

        iparser.close();

        return newTrainingData;
    }
}
