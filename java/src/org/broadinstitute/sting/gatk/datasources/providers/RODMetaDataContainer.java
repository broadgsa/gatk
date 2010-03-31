package org.broadinstitute.sting.gatk.datasources.providers;

import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
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
    private final HashMap<String, GATKFeature> nameMap = new HashMap<String, GATKFeature>();

    // we do allow duplicate class entries, so we need to store pairs of data
    private final List<Pair<Class, GATKFeature>> classMap = new ArrayList<Pair<Class, GATKFeature>>();

    public void addEntry(GATKFeature data) {
        nameMap.put(data.getName(),data);
        classMap.add(new Pair<Class, GATKFeature>(data.getClass(),data));
    }

    public Collection<GATKFeature> getSet(String name) {
        if (name == null) return nameMap.values();
        Set<GATKFeature> set = new HashSet<GATKFeature>();
        if (nameMap.containsKey(name)) set.add(nameMap.get(name));
        return set;
    }
    // the brute force (n) search ended up being faster than sorting and binary search in all but the most extreme cases (thousands of RODs at a location).
    public Collection<GATKFeature> getSet(Class cls) {
        Collection<GATKFeature> ret = new ArrayList<GATKFeature>();
        for (Pair<Class, GATKFeature> pair: classMap)
            if (pair.first.equals(cls)) ret.add(pair.second);
        return ret;
    }
}
