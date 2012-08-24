package org.broadinstitute.sting.utils.nanoScheduler;

/**
 * Created with IntelliJ IDEA.
 * User: depristo
 * Date: 8/24/12
 * Time: 9:57 AM
 * To change this template use File | Settings | File Templates.
 */
public class MapResult<MapType> implements Comparable<MapResult<MapType>> {
    final Integer id;
    final MapType value;

    public MapResult(final int id, final MapType value) {
        this.id = id;
        this.value = value;
    }

    public Integer getId() {
        return id;
    }

    public MapType getValue() {
        return value;
    }

    @Override
    public int compareTo(MapResult<MapType> o) {
        return getId().compareTo(o.getId());
    }
}
