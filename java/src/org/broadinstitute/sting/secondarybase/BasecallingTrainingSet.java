package org.broadinstitute.sting.secondarybase;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * BasecallingTrainingSet holds a set of raw read sequences, their raw intensities, and quality scores.
 *
 * @author Kiran Garimella
 */
public class BasecallingTrainingSet {
    private File bustardDir;
    private int lane;
    private int cycleBegin;
    private int cycleEnd;
    private int trainingLimit;

    private ArrayList<RawRead> trainingData;

    /**
     * Constructor for BasecallingTrainingSet.
     *
     * @param bustardDir     the Bustard directory for the sample
     * @param lane           the lane for the sample
     * @param cycleBegin     the start cycle for the beginning of the read (0-based, inclusive)
     * @param cycleEnd       the stop ccle for the end of the read (0-based, inclusive)
     * @param trainingLimit  the number of training reads to accept
     */
    public BasecallingTrainingSet(File bustardDir, int lane, int cycleBegin, int cycleEnd, int trainingLimit) {
        this.bustardDir = bustardDir;
        this.lane = lane;
        this.cycleBegin = cycleBegin;
        this.cycleEnd = cycleEnd;
        this.trainingLimit = trainingLimit;
    }

    /**
     * Get the training data array list
     *
     * @return  the arraylist of raw training reads
     */
    public ArrayList<RawRead> getTrainingData() {
        return this.trainingData;
    }

    /**
     * Set the training data array list
     *
     * @param trainingData  the arraylist of raw training reads
     */
    public void setTrainingData(ArrayList<RawRead> trainingData) {
        this.trainingData = trainingData;
    }

    /**
     * Take the first N reads that have no ambiguous bases and add them to the training set.
     */
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

    /**
     * Load a training set from perfect reads in an already-aligned bam file
     *
     * @param samIn      the SAM/BAM file to load the reads from
     * @param reference  the reference to which the reads should be compared
     */
    public void loadPreAlignedTrainingSet(File samIn, File reference) {
        Vector< HashMap<String, SAMRecord> > trainingReads = getPerfectAlignmentsByTile(samIn, reference);

        trainingData = correlateReadsAndIntensities(trainingReads);
    }

    /**
     * Find perfect reads and group them by tile.
     *
     * @param samIn      the SAM/BAM file to load the raeds from
     * @param reference  the reference to which the reads should be compared
     * @return a vector of perfect reads, grouped by tile
     */
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

    /**
     * Correlate the perfect reads with their raw intensities.  Sloooooooow.
     *
     * @param trainingReads  the perfect reads, grouped by tile
     * @return a training set of raw sequence, intensities, and quality scores (all set to 40 for these perfect bases)
     */
    private ArrayList<RawRead> correlateReadsAndIntensities(Vector<HashMap<String, SAMRecord>> trainingReads) {
        ArrayList<RawRead> newTrainingData = new ArrayList<RawRead>(trainingLimit);

        IlluminaParser iparser = new IlluminaParser(bustardDir, lane, cycleBegin, cycleEnd);

        int totalReadCount = 0;

        for (int tile = 1; tile < trainingReads.size(); tile++) {
            iparser.seekToTile(tile);

            int tileReadCount = 0;

            RawRead iread;
            while (trainingReads.get(tile) != null && tileReadCount < trainingReads.get(tile).size() && (iread = iparser.next()) != null) {
                String readKey = iread.getReadKey();

                if (trainingReads.get(tile).containsKey(readKey)) {
                    System.out.printf("\tTile %d: found %d of %d (%4.4f in tile, %4.4f total)            \r",
                                      tile,
                                      tileReadCount,
                                      trainingReads.get(tile).size(),
                                      ((double) tileReadCount)/((double) trainingReads.get(tile).size()),
                                      ((double) totalReadCount)/((double) trainingLimit));

                    byte[] quals = new byte[iread.getReadLength()];
                    for (int qualIndex = 0; qualIndex < quals.length; qualIndex++) {
                        quals[qualIndex] = 40;
                    }

                    iread.setQuals(quals);
                    newTrainingData.add(iread);

                    tileReadCount++;
                    totalReadCount++;
                }
            }
        }

        iparser.close();

        return newTrainingData;
    }
}
