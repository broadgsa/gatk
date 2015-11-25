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

package org.broadinstitute.gatk.utils.smithwaterman;

import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.collections.Pair;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Simple program to run SW performance test.
 *
 * // TODO -- should be replaced with Caliper before using again
 *
 * User: depristo
 * Date: 2/28/13
 * Time: 4:54 PM
 * To change this template use File | Settings | File Templates.
 */
public class SWPairwiseAlignmentMain {
    //    BELOW: main() method for testing; old implementations of the core methods are commented out below;
//           uncomment everything through the end of the file if benchmarking of new vs old implementations is needed.

    public static void main(String argv[]) {
//        String ref="CACGAGCATATGTGTACATGAATTTGTATTGCACATGTGTTTAATGCGAACACGTGTCATGTGTATGTGTTCACATGCATGTGTGTCT";
//        String read =   "GCATATGTTTACATGAATTTGTATTGCACATGTGTTTAATGCGAACACGTGTCATGTGTGTGTTCACATGCATGTG";

        String ref = null;
        String read = null;

        Map<String,List<String>> args = processArgs(argv);

        List<String> l = args.get("SEQ");
        args.remove("SEQ");
        if ( l == null ) {
            System.err.println("SEQ argument is missing. Two input sequences must be provided");
            System.exit(1);
        }
        if ( l.size() != 2 ) {
            System.err.println("Two input sequences (SEQ arguments) must be provided. Found "+l.size()+" instead");
            System.exit(1);
        }

        ref = l.get(0);
        read = l.get(1);

        Double m = extractSingleDoubleArg("MATCH",args);
        Double mm = extractSingleDoubleArg("MISMATCH",args);
        Double open = extractSingleDoubleArg("OPEN",args);
        Double ext = extractSingleDoubleArg("EXTEND",args);

        Boolean reverse = extractSingleBooleanArg("REVERSE",args);
        if ( reverse != null && reverse.booleanValue() == true ) {
            ref = Utils.reverse(ref);
            read = Utils.reverse(read);
        }

        Boolean print_mat = extractSingleBooleanArg("PRINT_MATRIX",args);
        Boolean cut = extractSingleBooleanArg("CUTOFF",args);
        if ( cut != null ) SWPairwiseAlignment.cutoff = cut;

        if ( args.size() != 0 ) {
            System.err.println("Unknown argument on the command line: "+args.keySet().iterator().next());
            System.exit(1);
        }

        final int w_match, w_mismatch, w_open, w_extend;

        w_match = (m == null ? 30 : m.intValue());
        w_mismatch = (mm == null ? -10 : mm.intValue());
        w_open = (open == null ? -10 : open.intValue());
        w_extend = (ext == null ? -2 : ext.intValue());


        SWPairwiseAlignment.keepScoringMatrix = true;
        SWPairwiseAlignment a = new SWPairwiseAlignment(ref.getBytes(),read.getBytes(),w_match,w_mismatch,w_open,w_extend);

        System.out.println("start="+a.getAlignmentStart2wrt1()+", cigar="+a.getCigar()+
                " length1="+ref.length()+" length2="+read.length());


        System.out.println();
        a.printAlignment(ref.getBytes(),read.getBytes());

        System.out.println();
        if ( print_mat != null && print_mat == true ) {
            print(a.SW,ref.getBytes(),read.getBytes());
        }
    }

    private static void print(final int[][] s, final byte[] a, final byte[] b) {
        int n = a.length+1;
        int m = b.length+1;
        System.out.print("         ");
        for ( int j = 1 ; j < m ; j++) System.out.printf(" %5c",(char)b[j-1]) ;
        System.out.println();

        for ( int i = 0 ; i < n ; i++) {
            if ( i > 0 ) System.out.print((char)a[i-1]);
            else System.out.print(' ');
            System.out.print("  ");
            for ( int j = 0; j < m ; j++ ) {
                System.out.printf(" %5.1f",s[i][j]);
            }
            System.out.println();
        }
    }


    static Pair<String,Integer> getArg(String prefix, String argv[], int i) {
        String arg = null;
        if ( argv[i].startsWith(prefix) ) {
            arg = argv[i].substring(prefix.length());
            if( arg.length() == 0 ) {
                i++;
                if ( i < argv.length ) arg = argv[i];
                else {
                    System.err.println("No value found after " + prefix + " argument tag");
                    System.exit(1);
                }
            }
            i++;
        }
        return new Pair<String,Integer>(arg,i);
    }

    static Map<String,List<String>> processArgs(String argv[]) {
        Map<String,List<String>> args = new HashMap<String,List<String>>();

        for ( int i = 0; i < argv.length ; i++ ) {
            String arg = argv[i];
            int pos = arg.indexOf('=');
            if ( pos < 0 ) {
                System.err.println("Argument "+arg+" is not of the form <ARG>=<VAL>");
                System.exit(1);
            }
            String val = arg.substring(pos+1);
            if ( val.length() == 0 ) {
                // there was a space between '=' and the value
                i++;
                if ( i < argv.length ) val = argv[i];
                else {
                    System.err.println("No value found after " + arg + " argument tag");
                    System.exit(1);
                }
            }
            arg = arg.substring(0,pos);

            List<String> l = args.get(arg);
            if ( l == null ) {
                l = new ArrayList<String>();
                args.put(arg,l);
            }
            l.add(val);
        }
        return args;
    }

    static Double extractSingleDoubleArg(String argname, Map<String,List<String>> args) {
        List<String> l = args.get(argname);
        args.remove(argname);
        if ( l == null ) return null;

        if ( l.size() > 1 ) {
            System.err.println("Only one "+argname+" argument is allowed");
            System.exit(1);
        }
        double d=0;
        try {
            d = Double.parseDouble(l.get(0));
        } catch ( NumberFormatException e) {
            System.err.println("Can not parse value provided for "+argname+" argument ("+l.get(0)+")");
            System.exit(1);
        }
        System.out.println("Argument "+argname+" set to "+d);
        return new Double(d);
    }


    static Boolean extractSingleBooleanArg(String argname, Map<String,List<String>> args) {
        List<String> l = args.get(argname);
        args.remove(argname);
        if ( l == null ) return null;

        if ( l.size() > 1 ) {
            System.err.println("Only one "+argname+" argument is allowed");
            System.exit(1);
        }
        if ( l.get(0).equals("true") ) return Boolean.valueOf(true);
        if ( l.get(0).equals("false") ) return Boolean.valueOf(false);
        System.err.println("Can not parse value provided for "+argname+" argument ("+l.get(0)+"); true/false are allowed");
        System.exit(1);
        return Boolean.valueOf(false); // This value isn't used because it is preceded by System.exit(1)
    }
}
