package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Oct 30, 2009
 */
public class NHashMap<T> extends HashMap<List<? extends Comparable<?>>, T> {

	private static final long serialVersionUID = 1L; //BUGBUG: what should I do here?

    

	public static <T extends Comparable<?>> List<T> makeList(T... args) {
        List<T> list = new ArrayList<T>();
        for (T arg : args)
        {
            list.add(arg);
        }
        return list;
    }
}





