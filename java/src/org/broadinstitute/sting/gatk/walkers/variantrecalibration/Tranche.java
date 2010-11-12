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
import org.broadinstitute.sting.utils.vcf.VCFUtils;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.*;
import java.util.*;

/**
 */

public class Tranche {
    public double fdr, pCut, knownTiTv, novelTiTv;
    public int numKnown,numNovel;
    public String name;

    public Tranche(double fdr, double pCut, int numKnown, double knownTiTv, int numNovel, double novelTiTv) {
        this(fdr, pCut, numKnown, knownTiTv, numNovel, novelTiTv, "anonymous");
    }

    public Tranche(double fdr, double pCut, int numKnown, double knownTiTv, int numNovel, double novelTiTv, String name) {
        this.fdr = fdr;
        this.pCut = pCut;
        this.novelTiTv = novelTiTv;
        this.numNovel = numNovel;
        this.knownTiTv = knownTiTv;
        this.numKnown = numKnown;
        this.name = name;
    }

    public String toString() {
        return String.format("[Tranche cut = %.3f with %d novels @ %.2f]", pCut, numNovel, novelTiTv);
    }

    public static String tranchesString(List<Tranche> tranches) {
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        PrintStream stream = new PrintStream(bytes);

        stream.println("FDRtranche,numKnown,numNovel,knownTiTv,novelTiTv,pCut,filterName");

        Tranche prev = null;
        for ( Tranche t : tranches ) {
            stream.printf("%.2f,%d,%d,%.4f,%.4f,%.4f,FDRtranche%.2fto%.2f%n",
                    t.fdr,t.numKnown,t.numNovel,t.knownTiTv,t.novelTiTv, t.pCut,
                    (prev == null ? 0.0 : prev.fdr), t.fdr);
            prev = t;
        }

        return bytes.toString();
    }

    private static double getDouble(Map<String,String> bindings, String key) {
        return bindings.containsKey(key) ? Double.valueOf(bindings.get(key)) : -1.0;
    }

    private static int getInteger(Map<String,String> bindings, String key) {
        return bindings.containsKey(key) ? Integer.valueOf(bindings.get(key)) : -1;
    }

    public static List<Tranche> readTraches(File f) {
        String[] header = null;
        List<Tranche> tranches = new ArrayList<Tranche>();

        try {
            for( final String line : new XReadLines(f) ) {
                final String[] vals = line.split(",");
                if( header == null ) {
                    header = vals;
                } else {
                    Map<String,String> bindings = new HashMap<String, String>();
                    for ( int i = 0; i < vals.length; i++ ) bindings.put(header[i], vals[i]);
                    tranches.add(new Tranche(getDouble(bindings,"FDRtranche"),
                            getDouble(bindings,"pCut"),
                            getInteger(bindings,"numKnown"),
                            getDouble(bindings,"knownTiTv"),
                            getInteger(bindings,"numNovel"),
                            Math.max(getDouble(bindings,"novelTiTv"), getDouble(bindings,"novelTITV")),
                            bindings.get("filterName")));
                }
            }

            return tranches;
        } catch( FileNotFoundException e ) {
            throw new UserException.CouldNotCreateOutputFile(f, e);
        }
    }
}