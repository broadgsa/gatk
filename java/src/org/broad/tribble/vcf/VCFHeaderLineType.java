package org.broad.tribble.vcf;

/**
 * the type encodings we use for fields in VCF header lines
 */
public enum VCFHeaderLineType {
    Integer, Float, String, Character, Flag;

    public Object convert(String value, VCFCompoundHeaderLine.SupportedHeaderLineType hlt) {
        switch (this) {
            case Integer:
                return java.lang.Integer.valueOf(value); // the java.lang is needed since we use Integer as a enum name
            case Float:
                return java.lang.Float.valueOf(value);
            case String:
                return value;
            case Character:
                if (value.length()!= 0)
                    throw new IllegalStateException("INFO_TYPE." + this + " requires fields of length 1, what was provided was " + value);
                return value;
            case Flag:
                 if (hlt.allowFlagValues)
                    return value.equals("0") ? false : true;
            default:
                throw new IllegalStateException("INFO_TYPE." + this + " doesn't have a set conversion approach");
        }
    }
}
