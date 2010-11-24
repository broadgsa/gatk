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

import org.broad.tribble.util.variantcontext.VariantContext;
import java.util.*;

/* Some methods for extracting RefSeq-related data from annotated VCF INFO fields:
 */
public class RefSeqDataParser {
    private static String REFSEQ_PREFIX = "refseq.";

    private static String NUM_RECORDS_KEY = REFSEQ_PREFIX + "numMatchingRecords";
    private static String NAME_KEY = REFSEQ_PREFIX + "name";
    private static String NAME2_KEY = REFSEQ_PREFIX + "name2";

    private static Map<String, String> getRefSeqEntriesToNames(VariantContext vc, boolean getName2) {
        String nameKeyToUse = getName2 ? NAME2_KEY : NAME_KEY;
        String nameKeyToUseMultiplePrefix = nameKeyToUse + "_";

        Map<String, String> entriesToNames = new HashMap<String, String>();
        Integer numRecords = vc.getAttributeAsIntegerNoException(NUM_RECORDS_KEY);
        if (numRecords != null) {
            for (int i = 1; i <= numRecords; i++) {
                String key = nameKeyToUseMultiplePrefix + i;
                String name = vc.getAttributeAsStringNoException(key);
                if (name != null)
                    entriesToNames.put(key, name);
            }
        }
        else { // no entry with the # of records:
            String name = vc.getAttributeAsStringNoException(nameKeyToUse);
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