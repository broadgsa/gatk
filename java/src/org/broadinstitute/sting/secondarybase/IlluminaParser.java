package org.broadinstitute.sting.secondarybase;

import edu.mit.broad.picard.util.BasicTextFileParser;
import edu.mit.broad.picard.util.PasteParser;
import org.broadinstitute.sting.utils.StingException;

import java.io.Closeable;
import java.io.File;
import java.io.FilenameFilter;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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

    public IlluminaParser(File bustardDir, int lane) {
        this.bustardDir = bustardDir;
        this.firecrestDir = bustardDir.getParentFile();
        this.lane = lane;

        initializeParser();
    }

    public IlluminaParser(File bustardDir, File firecrestDir, int lane) {
        this.bustardDir = bustardDir;
        this.firecrestDir = firecrestDir;
        this.lane = lane;

        initializeParser();
    }

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

    private FilenameFilter getFilenameFilter(final String suffix) {
        return new FilenameFilter() {
            public boolean accept(File file, String s) {
                Pattern pseq = Pattern.compile(String.format("s_%d_\\d+_%s\\.txt(?!.+old.+)(\\.gz)?", lane, suffix));
                Matcher mseq = pseq.matcher(s);

                return mseq.find();
            }
        };
    }

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

    public int numTiles() { return intfiles.length; }

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

    public boolean hasNext() {
        return (currentTileParser.hasNext() || seekToTile(currentTileIndex + 1));
    }

    public boolean next() {
        if (hasNext()) {
            currentParseResult = currentTileParser.next();

            return true;
        }

        return false;
    }

    public String[][] getCurrentParseResult() {
        return currentParseResult;
    }

    public void remove() {
        throw new UnsupportedOperationException("IlluminaParser.remove() method is not supported.");
    }

    public void close() {
        currentTileParser.close();
    }

    public RawRead getRawRead() {
        return getSubset(0, currentParseResult[1][4].length() - 1);
    }

    public RawRead getSubset(int cycleStart, int cycleStop) {
        return new RawRead(currentParseResult, cycleStart, cycleStop);
    }
}
