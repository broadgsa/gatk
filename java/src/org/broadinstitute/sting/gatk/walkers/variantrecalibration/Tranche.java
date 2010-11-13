/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.variantrecalibration;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.vcf.VCFUtils;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.*;
import java.util.*;

/**
 */

public class Tranche implements Comparable<Tranche> {
    private static final int CURRENT_VERSION = 2;

    public double fdr, pCut, targetTiTv, knownTiTv, novelTiTv;
    public int numKnown,numNovel;
    public String name;

    public Tranche(double fdr, double targetTiTv, double pCut, int numKnown, double knownTiTv, int numNovel, double novelTiTv) {
        this(fdr, targetTiTv, pCut, numKnown, knownTiTv, numNovel, novelTiTv, "anonymous");
    }

    public Tranche(double fdr, double targetTiTv, double pCut, int numKnown, double knownTiTv, int numNovel, double novelTiTv, String name) {
        this.fdr = fdr;
        this.targetTiTv = targetTiTv;
        this.pCut = pCut;
        this.novelTiTv = novelTiTv;
        this.numNovel = numNovel;
        this.knownTiTv = knownTiTv;
        this.numKnown = numKnown;
        this.name = name;

//        if ( targetTiTv < 0.5 || targetTiTv > 10 )
//            throw new UserException("Target Ti/Tv ratio is unreasonable " + targetTiTv);
//
//        if ( numKnown < 0 || numNovel < 0)
//            throw new ReviewedStingException("Invalid tranch - no. variants is < 0 : known " + numKnown + " novel " + numNovel);

        if ( name == null )
            throw new ReviewedStingException("BUG -- name cannot be null");

    }

    public int compareTo(Tranche other) {
        return Double.compare(this.fdr,  other.fdr);
    }

    public String toString() {
        return String.format("Tranche fdr=%.2f minQual=%.2f known=(%d @ %.2f) novel=(%d @ %.2f) name=%s]",
                fdr, pCut, numKnown, knownTiTv, numNovel, novelTiTv, name);
    }

    /**
     * Returns an appropriately formated string representing the raw tranches file on disk.
     *
     * @param rawTranches
     * @return
     */
    public static String tranchesString(List<Tranche> rawTranches) {
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        PrintStream stream = new PrintStream(bytes);

        List<Tranche> tranches = new ArrayList<Tranche>();
        tranches.addAll(rawTranches);
        Collections.sort(tranches);

        stream.println("# Variant quality score tranches file");
        stream.println("# Version number " + CURRENT_VERSION);
        stream.println("FDRtranche,targetTiTv,numKnown,numNovel,knownTiTv,novelTiTv,pCut,filterName");

        Tranche prev = null;
        for ( Tranche t : tranches ) {
            stream.printf("%.2f,%.2f,%d,%d,%.4f,%.4f,%.2f,FDRtranche%.2fto%.2f%n",
                    t.fdr,t.targetTiTv,t.numKnown,t.numNovel,t.knownTiTv,t.novelTiTv, t.pCut,
                    (prev == null ? 0.0 : prev.fdr), t.fdr);
            prev = t;
        }

        return bytes.toString();
    }

    private static double getDouble(Map<String,String> bindings, String key, boolean required) {
        if ( bindings.containsKey(key) )
            return Double.valueOf(bindings.get(key));
        else if ( required ) {
            throw new UserException.MalformedFile("Malformed tranches file.  Missing required key " + key);
        }
        else
            return -1;
    }

    private static int getInteger(Map<String,String> bindings, String key, boolean required) {
        if ( bindings.containsKey(key) )
            return Integer.valueOf(bindings.get(key));
        else if ( required ) {
            throw new UserException.MalformedFile("Malformed tranches file.  Missing required key " + key);
        }
        else
            return -1;
    }

    /**
     * Returns a list of tranches, sorted from most to least specific, read in from file f
     * 
     * @param f
     * @return
     */
    public static List<Tranche> readTraches(File f) {
        String[] header = null;
        List<Tranche> tranches = new ArrayList<Tranche>();

        try {
            for( final String line : new XReadLines(f) ) {
                if ( line.startsWith("#") )
                    continue;

                final String[] vals = line.split(",");
                if( header == null ) {
                    header = vals;
                } else {
                    if ( header.length != vals.length )
                        throw new UserException.MalformedFile(f, "Line had too few/many fields.  Header = " + header.length + " vals " + vals.length + " line " + line);

                    Map<String,String> bindings = new HashMap<String, String>();
                    for ( int i = 0; i < vals.length; i++ ) bindings.put(header[i], vals[i]);
                    tranches.add(new Tranche(getDouble(bindings,"FDRtranche", true),
                            getDouble(bindings,"targetTiTv", false),
                            getDouble(bindings,"pCut", true),
                            getInteger(bindings,"numKnown", false),
                            getDouble(bindings,"knownTiTv", false),
                            getInteger(bindings,"numNovel", true),
                            Math.max(getDouble(bindings,"novelTiTv", false), getDouble(bindings,"novelTITV", false)),
                            bindings.get("filterName")));
                }
            }

            Collections.sort(tranches);   // sort this in the standard order
            return tranches;
        } catch( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(f, e);
        }
    }
}
