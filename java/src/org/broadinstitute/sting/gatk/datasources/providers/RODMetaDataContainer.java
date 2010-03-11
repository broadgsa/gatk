package org.broadinstitute.sting.gatk.datasources.providers;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.utils.Pair;

import java.util.*;


/**
 * 
 * @author aaron 
 * 
 * Class RODMetaDataContainer
 *
 * stores both the name and the class for each ROD.  This class assumes that:
 *
 * -Names must be unique
 * -Classes are allowed to have dupplicates
 *
 * This class encapsulates the ref data associations, and provides lookup by name and by
 * class type.
 *
 */
public class RODMetaDataContainer {
    // we only allow non-dupplicate ROD names, a HashMap is fine
    private final HashMap<String, ReferenceOrderedDatum> nameMap = new HashMap<String, ReferenceOrderedDatum>();

    // we do allow duplicate class entries, so we need to store pairs of data
    private final List<Pair<Class, ReferenceOrderedDatum>> classMap = new ArrayList<Pair<Class, ReferenceOrderedDatum>>();

    public void addEntry(ReferenceOrderedDatum data) {
        nameMap.put(data.getName(),data);
        classMap.add(new Pair<Class, ReferenceOrderedDatum>(data.getClass(),data));
    }

    public Collection<ReferenceOrderedDatum> getSet(String name) {
        if (name == null) return nameMap.values();
        Set<ReferenceOrderedDatum> set = new HashSet<ReferenceOrderedDatum>();
        if (nameMap.containsKey(name)) set.add(nameMap.get(name));
        return set;
    }
    // the brute force (n) search ended up being faster than sorting and binary search in all but the most extreme cases (thousands of RODs at a location).
    public Collection<ReferenceOrderedDatum> getSet(Class cls) {
        Collection<ReferenceOrderedDatum> ret = new ArrayList<ReferenceOrderedDatum>();
        for (Pair<Class, ReferenceOrderedDatum> pair: classMap)
            if (pair.first.equals(cls)) ret.add(pair.second);
        return ret;
    }
}
