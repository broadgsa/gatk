package org.broadinstitute.sting.gatk.walkers.varianteval.util;

import java.util.TreeMap;

public class StateKey extends TreeMap<String, String> {
    public int hashCode() {
        int hashCode = 1;

        for (String key : this.keySet()) {
            String value = this.get(key);

            hashCode *= key.hashCode() + value.hashCode();
        }

        return hashCode;
    }

    public String toString() {
        String value = "";

        for ( String key : this.keySet() ) {
            //value += "\tstate " + key + ":" + this.get(key) + "\n";
            value += String.format("%s:%s;", key, this.get(key));
        }

        return value;
    }
}
