package org.broadinstitute.sting.gatk.dataSources.chunks;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: aaronmckenna
 * Date: Mar 29, 2009
 * Time: 8:35:16 PM
 * To change this template use File | Settings | File Templates.
 */
public class BasicDataShard<T> implements DataShard {

    List<T> list = new ArrayList<T>();
    int index = 0;

    public BasicDataShard(List<T> list) {
        this.list = list;
    }

    public boolean hasNext() {
        if (list.size() > index) {
            return true;
        }
        return false;
    }

    public T next() {
        return list.get(index);
    }

    public void remove() {
        list.remove(index);
    }
}
