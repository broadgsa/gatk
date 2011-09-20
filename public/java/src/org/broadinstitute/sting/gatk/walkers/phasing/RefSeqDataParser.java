/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

/* Some methods for extracting RefSeq-related data from annotated VCF INFO fields:
 */
public class RefSeqDataParser {
    private static String REFSEQ_PREFIX = "refseq.";

    private static String NUM_RECORDS_KEY = REFSEQ_PREFIX + "numMatchingRecords";
    private static String NAME_KEY = REFSEQ_PREFIX + "name";
    private static String NAME2_KEY = REFSEQ_PREFIX + "name2";

    private static String[] NAME_KEYS = {NAME_KEY, NAME2_KEY};

    private static Map<String, String> getRefSeqEntriesToNames(VariantContext vc, boolean getName2) {
        String nameKeyToUse = getName2 ? NAME2_KEY : NAME_KEY;
        String nameKeyToUseMultiplePrefix = nameKeyToUse + "_";

        Map<String, String> entriesToNames = new HashMap<String, String>();
        int numRecords = vc.getAttributeAsInt(NUM_RECORDS_KEY, -1);
        if (numRecords != -1) {
            boolean done = false;

            if (numRecords == 1) { // Check if perhaps the single record doesn't end with "_1":
                String name = vc.getAttributeAsString(nameKeyToUse, null);
                if (name != null) {
                    entriesToNames.put(nameKeyToUse, name);
                    done = true;
                }
            }

            if (!done) {
                for (int i = 1; i <= numRecords; i++) {
                    String key = nameKeyToUseMultiplePrefix + i;
                    String name = vc.getAttributeAsString(key, null);
                    if (name != null)
                        entriesToNames.put(key, name);
                }
            }
        }
        else { // no entry with the # of records:
            String name = vc.getAttributeAsString(nameKeyToUse, null);
            if (name != null) {
                entriesToNames.put(nameKeyToUse, name);
            }
            else { // Check all INFO fields for a match (if there are multiple entries):
                for (Map.Entry<String, Object> entry : vc.getAttributes().entrySet()) {
                    String key = entry.getKey();
                    if (key.startsWith(nameKeyToUseMultiplePrefix))
                        entriesToNames.put(key, entry.getValue().toString());
                }
            }
        }
        return entriesToNames;
    }

    private static Map<String, String> getRefSeqEntriesToNames(VariantContext vc) {
        return getRefSeqEntriesToNames(vc, false);
    }

    public static Set<String> getRefSeqNames(VariantContext vc, boolean getName2) {
        return new TreeSet<String>(getRefSeqEntriesToNames(vc, getName2).values());
    }

    public static Set<String> getRefSeqNames(VariantContext vc) {
        return getRefSeqNames(vc, false);
    }

    public static Map<String, Object> getMergedRefSeqNameAttributes(VariantContext vc1, VariantContext vc2) {
        Map<String, Object> refSeqNameAttribs = new HashMap<String, Object>();

        Map<String, RefSeqEntry> entriesMap1 = getAllRefSeqEntriesByName(vc1);
        Map<String, RefSeqEntry> entriesMap2 = getAllRefSeqEntriesByName(vc2);

        Set<String> commonNames = entriesMap1.keySet();
        commonNames.retainAll(entriesMap2.keySet());
        boolean addSuffix = commonNames.size() > 1;
        int nextCount = 1;

        for (String name : commonNames) {
            RefSeqEntry refseq1 = entriesMap1.get(name);
            RefSeqEntry refseq2 = entriesMap2.get(name);

            String keySuffix = "";
            if (addSuffix)
                keySuffix = "_" + nextCount;

            boolean added = false;
            for (String key : NAME_KEYS) {
                Object obj1 = refseq1.info.get(key);
                Object obj2 = refseq2.info.get(key);
                if (obj1 != null && obj2 != null && obj1.equals(obj2)) {
                    added = true;
                    String useKey = key + keySuffix;
                    refSeqNameAttribs.put(useKey, obj1);
                }
            }
            if (added)
                nextCount++;
        }
        int totalCount = nextCount - 1; // since incremented count one extra time
        if (totalCount > 1)
            refSeqNameAttribs.put(NUM_RECORDS_KEY, totalCount);

        return refSeqNameAttribs;
    }

    public static Map<String, Object> removeRefSeqAttributes(Map<String, Object> attributes) {
        Map<String, Object> removedRefSeqAttributes = new HashMap<String, Object>(attributes);

        Iterator<Map.Entry<String, Object>> attrIt = removedRefSeqAttributes.entrySet().iterator();
        while (attrIt.hasNext()) {
            String key = attrIt.next().getKey();
            if (key.startsWith(REFSEQ_PREFIX))
                attrIt.remove();
        }

        return removedRefSeqAttributes;
    }

    private static Map<String, RefSeqEntry> getAllRefSeqEntriesByName(VariantContext vc) {
        Map<String, RefSeqEntry> nameToEntries = new TreeMap<String, RefSeqEntry>();

        List<RefSeqEntry> allEntries = getAllRefSeqEntries(vc);
        for (RefSeqEntry entry : allEntries) {
            Object name = entry.info.get(NAME_KEY);
            if (name != null)
                nameToEntries.put(name.toString(), entry);
        }

        return nameToEntries;
    }

    // Returns a List of SEPARATE Map<refseq.ENTRY, refseq.VALUE> for EACH RefSeq annotation (i.e., each gene), stripping out the "_1", "_2", etc.
    private static List<RefSeqEntry> getAllRefSeqEntries(VariantContext vc) {
        List<RefSeqEntry> allRefSeq = new LinkedList<RefSeqEntry>();

        for (Map.Entry<String, String> entryToName : getRefSeqEntriesToNames(vc).entrySet()) {
            String entry = entryToName.getKey();
            String entrySuffix = entry.replaceFirst(NAME_KEY, "");
            allRefSeq.add(new RefSeqEntry(vc, entrySuffix));
        }

        return allRefSeq;
    }

    private static class RefSeqEntry {
        public Map<String, Object> info;

        public RefSeqEntry(VariantContext vc, String entrySuffix) {
            this.info = new HashMap<String, Object>();

            for (Map.Entry<String, Object> attribEntry : vc.getAttributes().entrySet()) {
                String key = attribEntry.getKey();
                if (key.startsWith(REFSEQ_PREFIX) && key.endsWith(entrySuffix)) {
                    String genericKey = key.replaceAll(entrySuffix, "");
                    this.info.put(genericKey, attribEntry.getValue());
                }
            }
        }
    }
}