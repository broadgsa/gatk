package org.broadinstitute.sting.secondarybase;

import edu.mit.broad.picard.util.BasicTextFileParser;
import edu.mit.broad.picard.util.PasteParser;
import org.broadinstitute.sting.utils.StingException;

import java.io.Closeable;
import java.io.File;
import java.io.FilenameFilter;
import java.util.Arrays;
import java.util.Comparator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * IlluminaParser parses raw Illumina data (raw intensities, basecalled read sequences and quality scores)
 * and presents it to the developer in an easy-to-use form.  It also permits random tile jumping.
 *
 * WARNING: This parser does not understand newer GAPipeline data formats, and instead relies on the older
 * data formats that may have been generated with older Illumina tools.  As a result, this may return
 * suboptimal data.  This parser only exists temporarily until the Picard team writes a much more sensible
 * version.  Proceed with caution.
 *
 * @author Kiran Garimella
 */
public class IlluminaParser implements Closeable {
    private File bustardDir;
    private File firecrestDir;
    private int lane;

    private File[] intfiles;
    private File[] seqfiles;
    private File[] prbfiles;

    private int currentTileIndex;
    private PasteParser currentTileParser;
    private String[][] currentParseResult;

    /**
     * Construct an IlluminaParser given the Bustard directory and lane.  Infer the Firecrest directory.
     *
     * @param bustardDir  the Illumina Bustard directory
     * @param lane        the Illumina lane
     */
    public IlluminaParser(File bustardDir, int lane) {
        this.bustardDir = bustardDir;
        this.firecrestDir = bustardDir.getParentFile();
        this.lane = lane;

        initializeParser();
    }

    /**
     * Construct an IlluminaParser given the Bustard directory, Firecrest directory and lane.
     *
     * @param bustardDir    the Illumina Bustard directory
     * @param firecrestDir  the Illumina Firecrest directory
     * @param lane          the Illumina lane
     */
    public IlluminaParser(File bustardDir, File firecrestDir, int lane) {
        this.bustardDir = bustardDir;
        this.firecrestDir = firecrestDir;
        this.lane = lane;

        initializeParser();
    }

    /**
     * Initialize the parser and seek to the first tile.
     */
    private void initializeParser() {
        intfiles = firecrestDir.listFiles(getFilenameFilter("int"));
        seqfiles = bustardDir.listFiles(getFilenameFilter("seq"));
        prbfiles = bustardDir.listFiles(getFilenameFilter("prb"));

        if (intfiles.length != seqfiles.length || intfiles.length != prbfiles.length || seqfiles.length != prbfiles.length) {
            throw new StingException(
                String.format("File list lengths are unequal (int:%d, seq:%d, prb:%d)",
                              intfiles.length,
                              seqfiles.length,
                              prbfiles.length)
            );
        }

        Arrays.sort(intfiles, getTileSortingComparator());
        Arrays.sort(seqfiles, getTileSortingComparator());
        Arrays.sort(prbfiles, getTileSortingComparator());

        seekToTile(1);

        // Todo: put some more consistency checks here
    }

    /**
     * Get the filename filter for files of a given type.
     *
     * @param suffix  the type (i.e. 'int', 'seq', 'prb').
     * @return the filename filter
     */
    private FilenameFilter getFilenameFilter(final String suffix) {
        return new FilenameFilter() {
            public boolean accept(File file, String s) {
                Pattern pseq = Pattern.compile(String.format("s_%d_\\d+_%s\\.txt(?!.+old.+)(\\.gz)?", lane, suffix));
                Matcher mseq = pseq.matcher(s);

                return mseq.find();
            }
        };
    }

    /**
     * Get a comparator that sorts by tile.
     *
     * @return the comparator that sorts by tile.
     */
    private Comparator<File> getTileSortingComparator() {
        return new Comparator<File>() {
            public int compare(File file1, File file2) {
                Pattern ptile = Pattern.compile(String.format("s_%d_(\\d+)_", lane));

                Matcher mtile1 = ptile.matcher(file1.getName());
                Matcher mtile2 = ptile.matcher(file2.getName());

                if (mtile1.find() && mtile2.find()) {
                    int tile1 = Integer.valueOf(mtile1.group(1));
                    int tile2 = Integer.valueOf(mtile2.group(1));

                    if (tile1 < tile2) { return -1; }
                    else if (tile1 > tile2) { return 1; }

                    return 0;
                }

                throw new StingException("Tile filenames ('" + file1.getName() + "' or '" + file2.getName() + "') did not match against regexp pattern ('" + ptile.pattern() + "')");
            }
        };
    }

    /**
     * Return the number of tiles.
     *
     * @return  the number of tiles.
     */
    public int numTiles() { return intfiles.length; }

    /**
     * Seek to a specified tile.
     *
     * @param tile  the tile to which we should seek
     * @return true if we were able to seek to the tile, false if otherwise
     */
    public boolean seekToTile(int tile) {
        if (tile < intfiles.length - 1) {
            currentTileIndex = tile - 1;

            BasicTextFileParser intparser = new BasicTextFileParser(true, intfiles[currentTileIndex]);
            BasicTextFileParser seqparser = new BasicTextFileParser(true, seqfiles[currentTileIndex]);
            BasicTextFileParser prbparser = new BasicTextFileParser(true, prbfiles[currentTileIndex]);

            currentTileParser = new PasteParser(intparser, seqparser, prbparser);

            return true;
        }

        return false;
    }

    /**
     * Returns whether the parser has any more data to go through.
     *
     * @return  true if there's data left, false if otherwise
     */
    public boolean hasNext() {
        return (currentTileParser.hasNext() || seekToTile(currentTileIndex + 1));
    }

    /**
     * Advance the parser to the next read.
     *
     * @return true if successful, false if otherwise
     */
    public boolean next() {
        if (hasNext()) {
            currentParseResult = currentTileParser.next();

            return true;
        }

        return false;
    }

    /**
     * Returns the result from the current parse as an matrix of Strings.
     *
     * @return the matrix of Strings containing the current parse result
     */
    public String[][] getCurrentParseResult() {
        return currentParseResult;
    }

    /**
     * Removes, um, something, but in reality, does nothing.
     */
    public void remove() {
        throw new UnsupportedOperationException("IlluminaParser.remove() method is not supported.");
    }

    /**
     * Close the current tile.
     */
    public void close() {
        currentTileParser.close();
    }

    /**
     * Returns a raw read containing the raw intensities, read sequence, and quality scores.
     *
     * @return the raw read
     */
    public RawRead getRawRead() {
        return getSubset(0, currentParseResult[1][4].length() - 1);
    }

    /**
     * Returns a subset of the current parse result as a raw read.
     *
     * @param cycleStart  the starting cycle for the desired subset
     * @param cycleStop   the ending cycle for the desired subset
     * @return the subset of the current parse result as a raw read
     */
    public RawRead getSubset(int cycleStart, int cycleStop) {
        return new RawRead(currentParseResult, cycleStart, cycleStop);
    }
}
