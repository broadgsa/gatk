package org.broadinstitute.sting.playground.piecemealannotator;

import edu.mit.broad.picard.util.BasicTextFileParser;
import edu.mit.broad.picard.util.PasteParser;
import org.broadinstitute.sting.secondarybase.RawRead;
import org.broadinstitute.sting.utils.StingException;

import java.io.Closeable;
import java.io.File;
import java.util.Iterator;

public class IlluminaTile implements Iterator<RawRead>, Iterable<RawRead>, Closeable {
    private PasteParser parser;

    private RawRead next;
    private boolean isIterating = false;

    public IlluminaTile(File bustardDir, int lane, int tile) {
        //BasicTextFileParser intparser = new BasicTextFileParser(true, new File(bustardDir.getParent() + "/" + String.format("s_%d_%04d_int.txt.gz", lane, tile)));
        //BasicTextFileParser seqparser = new BasicTextFileParser(true, new File(bustardDir.getAbsolutePath() + "/" + String.format("s_%d_%04d_seq.txt.gz", lane, tile)));
        //BasicTextFileParser prbparser = new BasicTextFileParser(true, new File(bustardDir.getAbsolutePath() + "/" + String.format("s_%d_%04d_prb.txt.gz", lane, tile)));
        
        BasicTextFileParser intparser = new BasicTextFileParser(true, findFile(bustardDir.getParentFile(), lane, tile, "int"));
        BasicTextFileParser seqparser = new BasicTextFileParser(true, findFile(bustardDir, lane, tile, "seq"));
        BasicTextFileParser prbparser = new BasicTextFileParser(true, findFile(bustardDir, lane, tile, "prb"));

        parser = new PasteParser(intparser, seqparser, prbparser);
    }

    private File findFile(File dir, int lane, int tile, String ext) {
        File file1 = new File(dir.getAbsolutePath() + "/" + String.format("s_%d_%04d_%s.txt.gz", lane, tile, ext));
        if (file1.exists()) { return file1; }

        File file2 = new File(dir.getAbsolutePath() + "/" + String.format("s_%d_%04d_%s.txt", lane, tile, ext));
        if (file2.exists()) { return file2; }

        throw new StingException(String.format("Can't find file '%s' or '%s'", file1.getName(), file2.getName()));

    }

    public boolean hasNext() {
        if (!isIterating) {
            iterator();
        }

        return parser.hasNext();
    }

    public RawRead next() {
        if (hasNext()) {
            next = new RawRead(parser.next());
            return next;
        }

        return null;
    }

    public void remove() {
        throw new StingException("Remove is not implemented by IlluminaTile");
    }

    public void close() {
        parser.close();
        isIterating = false;
    }

    public Iterator<RawRead> iterator() {
        if (isIterating) {
            throw new StingException("IlluminaTile is already iterating");
        }

        isIterating = true;
        next = new RawRead(parser.next());

        return this;
    }
}
