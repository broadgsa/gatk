package org.broadinstitute.sting.gatk.walkers.varianteval.util;

import java.util.Map;
import java.util.TreeMap;

public class StateKey extends TreeMap<String, String> {
//    public int hashCode() {
//        int hashCode = 1;
//
//        for (final Map.Entry<String,String> pair : this.entrySet()) {
//            hashCode *= pair.getKey().hashCode() + pair.getValue().hashCode();
//        }
//
//        return hashCode;
//    }

    public String toString() {
        String value = "";

        for ( final String key : this.keySet() ) {
            //value += "\tstate " + key + ":" + this.get(key) + "\n";
            value += String.format("%s:%s;", key, this.get(key));
        }

        return value;
    }
}
