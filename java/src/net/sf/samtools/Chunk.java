package net.sf.samtools;

import net.sf.picard.PicardException;

import java.util.List;
import java.util.ArrayList;

/**
 * Represents a chunk stolen from the BAM file.  Originally a private static inner class within
 * BAMFileIndex; now breaking it out so that the sharding system can use it.
 *
 * @author mhanna
 * @version 0.1
 */
public class Chunk implements Cloneable,Comparable<Chunk> {

    private long mChunkStart;
    private long mChunkEnd;

    public Chunk(final long start, final long end) {
        mChunkStart = start;
        mChunkEnd = end;
    }

    protected Chunk clone() {
        return new Chunk(mChunkStart,mChunkEnd);
    }

    public long getChunkStart() {
        return mChunkStart;
    }

    public void setChunkStart(final long value) {
        mChunkStart = value;
    }

    public long getChunkEnd() {
        return mChunkEnd;
    }

    public void setChunkEnd(final long value) {
        mChunkEnd = value;
    }

    /**
     * The list of chunks is often represented as an array of
     * longs where every even-numbered index is a start coordinate
     * and every odd-numbered index is a stop coordinate.  Convert
     * from that format back to a list of chunks.
     * @param coordinateArray List of chunks to convert.
     * @return A list of chunks.
     */
    public static List<Chunk> toChunkList(long[] coordinateArray) {
        if(coordinateArray.length % 2 != 0)
            throw new PicardException("Data supplied does not appear to be in coordinate array format.");

        // TODO: possibly also check for monotonically increasing; this seems to be an implicit requirement of this format.
        List<Chunk> chunkList = new ArrayList<Chunk>();
        for(int i = 0; i < coordinateArray.length; i += 2)
            chunkList.add(new Chunk(coordinateArray[i],coordinateArray[i+1]));    

        return chunkList;
    }

    /**
     * The list of chunks is often represented as an array of
     * longs where every even-numbered index is a start coordinate
     * and every odd-numbered index is a stop coordinate.
     * @param chunks List of chunks to convert.
     * @return A long array of the format described above.
     */
    public static long[] toCoordinateArray(List<Chunk> chunks) {
        long[] coordinateArray = new long[chunks.size()*2];
        int position = 0;
        for(Chunk chunk: chunks) {
            coordinateArray[position++] = chunk.getChunkStart();
            coordinateArray[position++] = chunk.getChunkEnd();
        }
        return coordinateArray;
    }

    public int compareTo(final Chunk chunk) {
        int result = Long.signum(mChunkStart - chunk.mChunkStart);
        if (result == 0) {
            result = Long.signum(mChunkEnd - chunk.mChunkEnd);
        }
        return result;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final Chunk chunk = (Chunk) o;

        if (mChunkEnd != chunk.mChunkEnd) return false;
        if (mChunkStart != chunk.mChunkStart) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = (int) (mChunkStart ^ (mChunkStart >>> 32));
        result = 31 * result + (int) (mChunkEnd ^ (mChunkEnd >>> 32));
        return result;
    }

    @Override
    public String toString() {
        return String.format("%d:%d-%d:%d",mChunkStart >> 16,mChunkStart & 0xFFFF,mChunkEnd >> 16,mChunkEnd & 0xFFFF);
    }
}
