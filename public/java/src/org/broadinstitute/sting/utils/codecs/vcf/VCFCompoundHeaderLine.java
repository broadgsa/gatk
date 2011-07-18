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

package org.broadinstitute.sting.utils.codecs.vcf;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * a base class for compound header lines, which include info lines and format lines (so far)
 */
public abstract class VCFCompoundHeaderLine extends VCFHeaderLine implements VCFNamedHeaderLine {
    public enum SupportedHeaderLineType {
        INFO(true), FORMAT(false);

        public final boolean allowFlagValues;
        SupportedHeaderLineType(boolean flagValues) {
            allowFlagValues = flagValues;
        }
    }

    // the field types
    private String name;
    private int count = -1;
    private VCFHeaderLineCount countType;
    private String description;
    private VCFHeaderLineType type;

    // access methods
    public String getName() { return name; }
    public String getDescription() { return description; }
    public VCFHeaderLineType getType() { return type; }
    public VCFHeaderLineCount getCountType() { return countType; }
    public int getCount() {
        if ( countType != VCFHeaderLineCount.INTEGER )
            throw new ReviewedStingException("Asking for header line count when type is not an integer");
        return count;
    }

    // utility method
    public int getCount(int numAltAlleles) {
        int myCount;
        switch ( countType ) {
            case INTEGER: myCount = count; break;
            case UNBOUNDED: myCount = -1; break;
            case A: myCount = numAltAlleles; break;
            case G: myCount = ((numAltAlleles + 1) * (numAltAlleles + 2) / 2); break;
            default: throw new ReviewedStingException("Unknown count type: " + countType);
        }
        return myCount;
    }

    public void setNumberToUnbounded() {
        countType = VCFHeaderLineCount.UNBOUNDED;
        count = -1;
    }

    // our type of line, i.e. format, info, etc
    private final SupportedHeaderLineType lineType;

    /**
     * create a VCF format header line
     *
     * @param name         the name for this header line
     * @param count        the count for this header line
     * @param type         the type for this header line
     * @param description  the description for this header line
     * @param lineType     the header line type
     */
    protected VCFCompoundHeaderLine(String name, int count, VCFHeaderLineType type, String description, SupportedHeaderLineType lineType) {
        super(lineType.toString(), "");
        this.name = name;
        this.countType = VCFHeaderLineCount.INTEGER;
        this.count = count;
        this.type = type;
        this.description = description;
        this.lineType = lineType;
        validate();
    }

    /**
     * create a VCF format header line
     *
     * @param name         the name for this header line
     * @param count        the count type for this header line
     * @param type         the type for this header line
     * @param description  the description for this header line
     * @param lineType     the header line type
     */
    protected VCFCompoundHeaderLine(String name, VCFHeaderLineCount count, VCFHeaderLineType type, String description, SupportedHeaderLineType lineType) {
        super(lineType.toString(), "");
        this.name = name;
        this.countType = count;
        this.type = type;
        this.description = description;
        this.lineType = lineType;
        validate();
    }

    /**
     * create a VCF format header line
     *
     * @param line   the header line
     * @param version      the VCF header version
     * @param lineType     the header line type
     *
     */
    protected VCFCompoundHeaderLine(String line, VCFHeaderVersion version, SupportedHeaderLineType lineType) {
        super(lineType.toString(), "");
        Map<String,String> mapping = VCFHeaderLineTranslator.parseLine(version,line, Arrays.asList("ID","Number","Type","Description"));
        name = mapping.get("ID");
        count = -1;
        final String numberStr = mapping.get("Number");
        if ( numberStr.equals(VCFConstants.PER_ALLELE_COUNT) ) {
            countType = VCFHeaderLineCount.A;
        } else if ( numberStr.equals(VCFConstants.PER_GENOTYPE_COUNT) ) {
            countType = VCFHeaderLineCount.G;
        } else if ( ((version == VCFHeaderVersion.VCF4_0 || version == VCFHeaderVersion.VCF4_1) &&
                     numberStr.equals(VCFConstants.UNBOUNDED_ENCODING_v4)) ||
                    ((version == VCFHeaderVersion.VCF3_2 || version == VCFHeaderVersion.VCF3_3) &&
                     numberStr.equals(VCFConstants.UNBOUNDED_ENCODING_v3)) ) {
            countType = VCFHeaderLineCount.UNBOUNDED;
        } else {
            countType = VCFHeaderLineCount.INTEGER;
            count = Integer.valueOf(numberStr);

        }
        type = VCFHeaderLineType.valueOf(mapping.get("Type"));
        if (type == VCFHeaderLineType.Flag && !allowFlagValues())
            throw new IllegalArgumentException("Flag is an unsupported type for this kind of field");

        description = mapping.get("Description");
        if ( description == null && ALLOW_UNBOUND_DESCRIPTIONS ) // handle the case where there's no description provided
            description = UNBOUND_DESCRIPTION;
        
        this.lineType = lineType;

        validate();
    }

    private void validate() {
        if ( name == null || type == null || description == null || lineType == null )
            throw new IllegalArgumentException(String.format("Invalid VCFCompoundHeaderLine: key=%s name=%s type=%s desc=%s lineType=%s", 
                    super.getKey(), name, type, description, lineType ));
    }

    /**
     * make a string representation of this header line
     * @return a string representation
     */
    protected String toStringEncoding() {
        Map<String,Object> map = new LinkedHashMap<String,Object>();
        map.put("ID", name);
        Object number;
        switch ( countType ) {
            case A: number = VCFConstants.PER_ALLELE_COUNT; break;
            case G: number = VCFConstants.PER_GENOTYPE_COUNT; break;
            case UNBOUNDED: number = VCFConstants.UNBOUNDED_ENCODING_v4; break;
            case INTEGER:
            default: number = count;
        }
        map.put("Number", number);
        map.put("Type", type);
        map.put("Description", description);
        return lineType.toString() + "=" + VCFHeaderLine.toStringEncoding(map);
    }

    /**
     * returns true if we're equal to another compounder header line
     * @param o a compound header line
     * @return true if equal
     */
    public boolean equals(Object o) {
        if ( !(o instanceof VCFCompoundHeaderLine) )
            return false;
        VCFCompoundHeaderLine other = (VCFCompoundHeaderLine)o;
        return equalsExcludingDescription(other) &&
                description.equals(other.description);
    }

    public boolean equalsExcludingDescription(VCFCompoundHeaderLine other) {
        return count == other.count &&
                countType == other.countType &&
                type == other.type &&
                lineType == other.lineType &&
                name.equals(other.name);
    }

    public boolean sameLineTypeAndName(VCFCompoundHeaderLine other) {
        return  lineType == other.lineType &&
                name.equals(other.name);
    }

    /**
     * do we allow flag (boolean) values? (i.e. booleans where you don't have specify the value, AQ means AQ=true)
     * @return true if we do, false otherwise
     */
    abstract boolean allowFlagValues();

}
