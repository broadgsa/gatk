package net.sf.samtools;

/**
 * Represents a chunk stolen from the BAM file.  Originally a private static inner class within
 * BAMFileIndex; now breaking it out so that the sharding system can use it.
 *
 * @author mhanna
 * @version 0.1
 */
class Chunk implements Comparable<Chunk> {

    private long mChunkStart;
    private long mChunkEnd;

    Chunk(final long start, final long end) {
        mChunkStart = start;
        mChunkEnd = end;
    }

    long getChunkStart() {
        return mChunkStart;
    }

    void setChunkStart(final long value) {
        mChunkStart = value;
    }

    long getChunkEnd() {
        return mChunkEnd;
    }

    void setChunkEnd(final long value) {
        mChunkEnd = value;
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
}
