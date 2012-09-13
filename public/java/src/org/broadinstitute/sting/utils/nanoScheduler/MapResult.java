package org.broadinstitute.sting.utils.nanoScheduler;

/**
 * Holds the results of a map job suitable for producer/consumer threading
 * via a BlockingQueue
 */
class MapResult<MapType> extends EOFMarkedValue<MapType> implements Comparable<MapResult<MapType>> {
    final int jobID;

    /**
     * Create a new MapResult with value datum and jod jobID ID
     *
     * @param datum the value produced by the map job
     * @param jobID the id of the map job (for correctness testing)
     */
    MapResult(final MapType datum, final int jobID) {
        super(datum);
        this.jobID = jobID;
        if ( jobID < 0 ) throw new IllegalArgumentException("JobID must be >= 0");
    }

    MapResult(final int jobID) {
        super();
        this.jobID = jobID;
        if ( jobID < 0 ) throw new IllegalArgumentException("JobID must be >= 0");
    }

    /**
     * Create the EOF marker version of MapResult
     */
    MapResult() {
        super();
        this.jobID = Integer.MAX_VALUE;
    }

    /**
     * @return the job ID of the map job that produced this MapResult
     */
    public int getJobID() {
        return jobID;
    }

    /**
     * Compare these MapResults in order of JobID.
     *
     * @param o
     * @return
     */
    @Override
    public int compareTo(MapResult<MapType> o) {
        return Integer.valueOf(jobID).compareTo(o.getJobID());
    }

    @Override
    public String toString() {
        return "[MapResult id=" + jobID + "]";
    }
}
