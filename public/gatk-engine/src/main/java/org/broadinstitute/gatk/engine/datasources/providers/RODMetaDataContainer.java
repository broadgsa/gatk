/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.datasources.providers;

import org.broadinstitute.gatk.utils.refdata.utils.GATKFeature;
import org.broadinstitute.gatk.utils.collections.Pair;

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
 * -Classes are allowed to have duplicates
 *
 * This class encapsulates the ref data associations, and provides lookup by name and by
 * class type.
 *
 */
public class RODMetaDataContainer {
    // we only allow non-duplicate ROD names, a HashMap is fine
    private final HashMap<String, GATKFeature> nameMap = new HashMap<String, GATKFeature>();

    // we do allow duplicate class entries, so we need to store pairs of data
    private final List<Pair<Class, GATKFeature>> classMap = new ArrayList<Pair<Class, GATKFeature>>();

    public void addEntry(GATKFeature data) {
        nameMap.put(data.getName(),data);
        classMap.add(new Pair<Class, GATKFeature>(data.getClass(),data));
    }

    public Collection<GATKFeature> getSet(String name) {
        if (name == null) return getSet();
        Set<GATKFeature> set = new HashSet<GATKFeature>();
        if (nameMap.containsKey(name)) set.add(nameMap.get(name));
        return set;
    }

    /**
     * get the feature contents of this container; the unfiltered set without their name association
     * @return
     */
    public Collection<GATKFeature> getSet() {
        return new ArrayList<GATKFeature>(nameMap.values());
    }

    // the brute force (n) search ended up being faster than sorting and binary search in all but the most extreme cases (thousands of RODs at a location).
    public Collection<GATKFeature> getSet(Class cls) {
        Collection<GATKFeature> ret = new ArrayList<GATKFeature>();
        for (Pair<Class, GATKFeature> pair: classMap)
            if (pair.first.equals(cls)) ret.add(pair.second);
        return ret;
    }
}
