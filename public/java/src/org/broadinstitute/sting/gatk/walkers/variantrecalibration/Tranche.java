/*
 * Copyright (c) 2011 The Broad Institute
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

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Mar 10, 2011
 */

public class Tranche implements Comparable<Tranche> {
    private static final int CURRENT_VERSION = 4;

    public double ts, minVQSLod, knownTiTv, novelTiTv;
    public int numKnown,numNovel;
    public String name;

    int accessibleTruthSites = 0;
    int callsAtTruthSites = 0;

    public Tranche(double ts, double minVQSLod, int numKnown, double knownTiTv, int numNovel, double novelTiTv, int accessibleTruthSites, int callsAtTruthSites) {
        this(ts, minVQSLod, numKnown, knownTiTv, numNovel, novelTiTv, accessibleTruthSites, callsAtTruthSites, "anonymous");
    }

    public Tranche(double ts, double minVQSLod, int numKnown, double knownTiTv, int numNovel, double novelTiTv, int accessibleTruthSites, int callsAtTruthSites, String name ) {
        this.ts = ts;
        this.minVQSLod = minVQSLod;
        this.novelTiTv = novelTiTv;
        this.numNovel = numNovel;
        this.knownTiTv = knownTiTv;
        this.numKnown = numKnown;
        this.name = name;

        this.accessibleTruthSites = accessibleTruthSites;
        this.callsAtTruthSites = callsAtTruthSites;

        if ( ts < 0.0 || ts > 100.0)
            throw new UserException("Target FDR is unreasonable " + ts);

        if ( numKnown < 0 || numNovel < 0)
            throw new ReviewedStingException("Invalid tranche - no. variants is < 0 : known " + numKnown + " novel " + numNovel);

        if ( name == null )
            throw new ReviewedStingException("BUG -- name cannot be null");
    }

    private double getTruthSensitivity() {
        return accessibleTruthSites > 0 ? callsAtTruthSites / (1.0*accessibleTruthSites) : 0.0;
    }

    public int compareTo(Tranche other) {
        return Double.compare(this.ts,  other.ts);
    }

    public String toString() {
        return String.format("Tranche ts=%.2f minVQSLod=%.4f known=(%d @ %.4f) novel=(%d @ %.4f) truthSites(%d accessible, %d called), name=%s]",
                ts, minVQSLod, numKnown, knownTiTv, numNovel, novelTiTv, accessibleTruthSites, callsAtTruthSites, name);
    }

    /**
     * Returns an appropriately formatted string representing the raw tranches file on disk.
     *
     * @param tranches
     * @return
     */
    public static String tranchesString( final List<Tranche> tranches ) {
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        final PrintStream stream = new PrintStream(bytes);

        Collections.sort(tranches);

        stream.println("# Variant quality score tranches file");
        stream.println("# Version number " + CURRENT_VERSION);
        stream.println("targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,accessibleTruthSites,callsAtTruthSites,truthSensitivity");

        Tranche prev = null;
        for ( Tranche t : tranches ) {
            stream.printf("%.2f,%d,%d,%.4f,%.4f,%.4f,TruthSensitivityTranche%.2fto%.2f,%d,%d,%.4f%n",
                    t.ts, t.numKnown, t.numNovel, t.knownTiTv, t.novelTiTv, t.minVQSLod,
                    (prev == null ? 0.0 : prev.ts), t.ts, t.accessibleTruthSites, t.callsAtTruthSites, t.getTruthSensitivity());
            prev = t;
        }

        return bytes.toString();
    }

    private static double getDouble(Map<String,String> bindings, String key, boolean required) {
        if ( bindings.containsKey(key) ) {
            String val = bindings.get(key);
            return Double.valueOf(val);
        }
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
    public static List<Tranche> readTranches(File f) {
        String[] header = null;
        List<Tranche> tranches = new ArrayList<Tranche>();

        try {
            for( final String line : new XReadLines(f) ) {
                if ( line.startsWith("#") )
                    continue;

                final String[] vals = line.split(",");
                if( header == null ) {
                    header = vals;
                    if ( header.length == 5 || header.length == 8 || header.length == 11 )
                        // old style tranches file, throw an error
                        throw new UserException.MalformedFile(f, "Unfortunately, your tranches file is from a previous version of this tool and cannot be used with the latest code.  Please rerun VariantRecalibrator");
                    if ( header.length != 10 )
                        throw new UserException.MalformedFile(f, "Expected 10 elements in header line " + line);
                } else {
                    if ( header.length != vals.length )
                        throw new UserException.MalformedFile(f, "Line had too few/many fields.  Header = " + header.length + " vals " + vals.length + ". The line was: " + line);

                    Map<String,String> bindings = new HashMap<String, String>();
                    for ( int i = 0; i < vals.length; i++ ) bindings.put(header[i], vals[i]);
                    tranches.add(new Tranche(getDouble(bindings,"targetTruthSensitivity", true),
                            getDouble(bindings,"minVQSLod", true),
                            getInteger(bindings,"numKnown", false),
                            getDouble(bindings,"knownTiTv", false),
                            getInteger(bindings,"numNovel", true),
                            getDouble(bindings,"novelTiTv", true),
                            getInteger(bindings,"accessibleTruthSites", false),
                            getInteger(bindings,"callsAtTruthSites", false),
                            bindings.get("filterName")));
                }
            }

            Collections.sort(tranches);
            return tranches;
        } catch( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(f, e);
        }
    }
}
