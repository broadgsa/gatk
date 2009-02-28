package edu.mit.broad.picard.directed;

import edu.mit.broad.picard.util.Interval;
import edu.mit.broad.picard.util.FormatUtil;
import edu.mit.broad.picard.io.IoUtil;
import edu.mit.broad.picard.PicardException;
import edu.mit.broad.sam.SAMFileHeader;
import edu.mit.broad.sam.SAMTextHeaderCodec;
import edu.mit.broad.sam.util.StringLineReader;

import java.util.*;
import java.io.*;

/**
 * Represents a list of intervals against a reference sequence that can be written to
 * and read from a file.  The file format is relatively simple and reflects the SAM
 * alignment format to a degree.
 *
 * A SAM style header must be present in the file which lists the sequence records
 * against which the intervals are described.  After the header the file then contains
 * records one per line in text format with the following values tab-separated:
 *   - Sequence name
 *   - Start position (1-based)
 *   - End position (1-based, end inclusive)
 *   - Strand (either + or -)
 *   - Interval name (an, ideally unique, name for the interval)
 *
 * @author Tim Fennell
 */
public class IntervalList implements Iterable<Interval> {
    private SAMFileHeader header;
    private List<Interval> intervals = new ArrayList<Interval>();

    /** Constructs a new interval list using the supplied header information. */
    public IntervalList(SAMFileHeader header) {
        if (header == null) {
            throw new IllegalArgumentException("SAMFileHeader must be supplied.");
        }
        this.header = header;
    }

    /** Gets the header (if there is one) for the interval list. */
    public SAMFileHeader getHeader() { return header; }

    /** Returns an iterator over the intervals. */
    public Iterator<Interval> iterator() { return this.intervals.iterator(); }

    /** Adds an interval to the list of intervals. */
    public void add(Interval interval) { this.intervals.add(interval); }

    /** Sorts the internal collection of intervals by coordinate. */
    public void sort() {
        Collections.sort(this.intervals, new IntervalCoordinateComparator(this.header));
        this.header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
    }

    /** Gets the set of intervals as held internally. */
    public List<Interval> getIntervals() {
        return Collections.unmodifiableList(this.intervals);
    }

    /**
     * Merges the list of intervals and then reduces them down where regions overlap
     * or are directly adjacent to one another.  During this process the "merged" interval
     * will retain the strand and name of the 5' most interval merged.
     *
     * @return the set of unique intervals condensed from the contained intervals
     */
    public List<Interval> getUniqueIntervals() {
        List<Interval> unique = new ArrayList<Interval>();
        ListIterator<Interval> iterator = this.intervals.listIterator();
        Interval previous = iterator.next();

        while (iterator.hasNext()) {
            Interval next = iterator.next();
            if (previous.intersects(next) || previous.abuts(next)) {
                previous = new Interval(previous.getSequence(),
                                        previous.getStart(),
                                        Math.max(previous.getEnd(), next.getEnd()),
                                        previous.isNegativeStrand(),
                                        previous.getName());
            }
            else {
                unique.add(previous);
                previous = next;
            }
        }

        if (previous != null) unique.add(previous);

        return unique;
    }

    /** Gets the (potentially redundant) sum of the length of the intervals in the list. */
    public long getBaseCount() {
        return Interval.countBases(this.intervals);
    }

    /** Gets the count of unique bases represented by the intervals in the list. */
    public long getUniqueBaseCount() {
        return Interval.countBases(getUniqueIntervals());
    }

    /**
     * Parses an interval list from a file.
     * @param file the file containing the intervals
     * @return an IntervalList object that contains the headers and intervals from the file
     */
    public static IntervalList fromFile(File file) {
        BufferedReader in = new BufferedReader(new InputStreamReader(IoUtil.openFileForReading(file)));

        try {
            // Setup a reader and parse the header
            StringBuilder builder = new StringBuilder(4096);
            String line = null;

            while ((line = in.readLine()) != null) {
                if (line.startsWith("@")) {
                    builder.append(line).append('\n');
                }
                else {
                    break;
                }
            }

            if (builder.length() == 0) {
                throw new IllegalStateException("Interval list file must contain header: " + file.getAbsolutePath());
            }

            StringLineReader headerReader = new StringLineReader(builder.toString());
            SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
            IntervalList list = new IntervalList(codec.decode(headerReader, file));
            
            // Then read in the intervals
            FormatUtil format = new FormatUtil();
            do {
                if (line.trim().length() == 0) continue; // skip over blank lines

                // Make sure we have the right number of fields
                String fields[] = line.split("\t");
                if (fields.length != 5) {
                    throw new PicardException("Invalid interval record contains " +
                                              fields.length + " fields: " + line);
                }

                // Then parse them out
                String seq = fields[0];
                int start = format.parseInt(fields[1]);
                int end   = format.parseInt(fields[2]);

                boolean negative;
                if (fields[3].equals("-")) negative = true;
                else if (fields[3].equals("+")) negative = false;
                else throw new IllegalArgumentException("Invalid strand field: " + fields[3]);

                String name = fields[4];

                Interval interval = new Interval(seq, start, end, negative, name);
                list.intervals.add(interval);
            }
            while ((line = in.readLine()) != null);

            return list;
        }
        catch (IOException ioe) {
            throw new PicardException("Error parsing interval list file: " + file.getAbsolutePath(), ioe);
        }
        finally {
            try { in.close(); } catch (Exception e) { /* do nothing */ }
        }
    }

    /**
     * Writes out the list of intervals to the supplied file.
     * @param file a file to write to.  If exists it will be overwritten.
     */
    public void write(File file) {
        try {
            BufferedWriter out = new BufferedWriter(new OutputStreamWriter(IoUtil.openFileForWriting(file)));
            FormatUtil format = new FormatUtil();

            // Write out the header
            if (this.header != null) {
                SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
                codec.encode(out, this.header);
            }

            // Write out the intervals
            for (Interval interval : this) {
                out.write(interval.getSequence());
                out.write('\t');
                out.write(format.format(interval.getStart()));
                out.write('\t');
                out.write(format.format(interval.getEnd()));
                out.write('\t');
                out.write(interval.isPositiveStrand() ? '+' : '-');
                out.write('\t');
                out.write(interval.getName());
                out.newLine();
            }

            out.flush();
            out.close();
        }
        catch (IOException ioe) {
            throw new PicardException("Error writing out interval list to file: " + file.getAbsolutePath(), ioe);
        }
    }
}

/**
 * Comparator that orders intervals based on their sequence index, by coordinate
 * then by strand and finally by name.
 */
class IntervalCoordinateComparator implements Comparator<Interval> {
    private SAMFileHeader header;

    /** Constructs a comparator using the supplied sequence header. */
    IntervalCoordinateComparator(SAMFileHeader header) {
        this.header = header;
    }

    public int compare(Interval lhs, Interval rhs) {
        int lhsIndex = this.header.getSequenceIndex(lhs.getSequence());
        int rhsIndex = this.header.getSequenceIndex(rhs.getSequence());
        int retval = lhsIndex - rhsIndex;

        if (retval == 0) retval = lhs.getStart() - rhs.getStart();
        if (retval == 0) retval = lhs.getEnd()   - rhs.getEnd();
        if (retval == 0) {
            if (lhs.isPositiveStrand() && rhs.isNegativeStrand()) retval = -1;
            else if (lhs.isNegativeStrand() && rhs.isPositiveStrand()) retval = 1;
        }
        if (retval == 0) {
            retval = lhs.getName().compareTo(rhs.getName());
        }

        return retval;
    }
}