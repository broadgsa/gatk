package edu.mit.broad.picard.directed;

import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Option;
import edu.mit.broad.picard.io.IoUtil;
import edu.mit.broad.picard.util.BasicTextFileParser;
import edu.mit.broad.picard.util.Interval;
import edu.mit.broad.picard.util.FormatUtil;
import edu.mit.broad.sam.SAMFileReader;
import edu.mit.broad.sam.SAMFileHeader;
import edu.mit.broad.sam.SAMSequenceRecord;

import java.io.File;
import java.util.List;

/**
 * Converts an arachne style map file to the new interval list format.
 *
 * @author Tim Fennell
 */
public class ArachneMapToIntervalList extends CommandLineProgram {
    @Option(shortName="M", doc="The path to an archne style map file") public File MAP;
    @Option(shortName="SD", doc="A sequence dictionary in SAM or BAM format") public File SEQUENCE_DICTIONARY;
    @Option(shortName="O", doc="The output file to write the interval list to") public File OUTPUT;
    @Option(shortName="P", doc="Prefix to use when generating names") public String PREFIX;

    /** Stock main method. */
    public static void main(String[] argv) {
        System.exit(new ArachneMapToIntervalList().instanceMain(argv));
    }

    protected int doWork() {
        IoUtil.assertFileIsReadable(MAP);
        IoUtil.assertFileIsReadable(SEQUENCE_DICTIONARY);
        IoUtil.assertFileIsWritable(OUTPUT);

        SAMFileReader sam = new SAMFileReader(SEQUENCE_DICTIONARY);
        SAMFileHeader header = sam.getFileHeader();
        List<SAMSequenceRecord> seqs = header.getSequences();
        IntervalList list = new IntervalList(header);

        BasicTextFileParser parser = new BasicTextFileParser(true, 3, MAP);
        FormatUtil format = new FormatUtil();
        int i=1;

        while (parser.hasNext()) {
            String[] fields = parser.next();
            int seqIndex = format.parseInt(fields[0]);
            int start    = format.parseInt(fields[1]) + 1;
            int end      = format.parseInt(fields[2]) + 1;
            String seq = seqs.get(seqIndex).getSequenceName();

            Interval interval = new Interval(seq, start, end, false, PREFIX + "_" + i++);
            list.add(interval);
        }

        list.sort();
        list.write(OUTPUT);

        return 0;
    }
}
