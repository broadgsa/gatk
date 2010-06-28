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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.qc;

import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.collections.PrimitivePair;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.commandline.Argument;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMReadGroupRecord;

import java.util.*;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedWriter;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Apr 9, 2010
 * Time: 12:16:41 PM
 * To change this template use File | Settings | File Templates.
 */

/**
 * Walks over the input data set, calculating the number of reads seen for diagnostic purposes.
 * Can also count the number of reads matching a given criterion using read filters (see the
 * --read-filter command line argument).  Simplest example of a read-backed analysis.
 */
@Requires({DataSource.READS})
public class CycleQualityWalker extends ReadWalker<Integer,Integer> {
    @Argument(fullName="mappedOnly", shortName="mo", doc="when this flag is set (default), statistics will be collected "+
                "on mapped reads only, while unmapped reads will be discarded", required=false)
    protected boolean MAPPED_ONLY = true;
    @Argument(fullName="maxReadLength", shortName="rl", doc="maximum read length", required=false)
    protected int MAX_READ_LENGTH = 500;
    @Argument(fullName="out_prefix",shortName="p",doc="prefix for output stat files",required=true)
    protected String PREFIX = null;
    @Argument(fullName="html",shortName="html",doc="produce html-formatted output (starting with h3-level tags) rather than plain text",required=false)
    protected boolean HTML = false;
    @Argument(fullName="qualThreshold", shortName="Q",doc="flag as problematic all cycles with av. qualities below the threshold",required=false)
    protected double QTHRESHOLD = 10.0;
    private Map<String,CycleStats[]> cyclesByLaneMap = null;
    private Map<String,CycleStats[]> cyclesByLibraryMap = null;

    public void initialize() {
        if ( PREFIX == null ) throw new StingException("Prefix for output file(s) must be specified");
        cyclesByLaneMap = new HashMap<String,CycleStats[]>();
        cyclesByLibraryMap = new HashMap<String,CycleStats[]>();
    }

    /** A trivial shortcut with the only purpose of avoiding copying the same code over and over again and
     * cluttering the main method. Takes lane (platform unit) and library strings already extracted from the read,
     * as well as read's base qualities and adds those qualities to the appropriate records of the specified by lane
     * and by library maps. The 'end' argument should be 1 or 2 for the first/second read in the pair, respectively.
     * If the read is unpaired (single end run), 'end' should be set to 0.
     * Being a shortcut, this method does not perform any checking except for adding a new
     * record to the respective map if it does not already contain a record for the specified lane or library.
     *
     * @param lane
     * @param library
     * @param quals
     */
    private void recordQuals(String lane, String library, byte [] quals, int end, Map<String,CycleStats[]> laneMap, Map<String,CycleStats[]> libMap) {

        CycleStats[] byLane = laneMap.get(lane);
        CycleStats[] byLib = libMap.get(library);

        // if end == 0 (single end lane), we allocate array of length 1, otherwise we need two
        // elements in the array in order to be able to collect statistics for each end in the pair independently
        if ( byLane == null ) laneMap.put(lane,byLane = new CycleStats[(end==0?1:2)]);
        if ( byLib == null ) libMap.put(library, byLib =new CycleStats[2]);

        if ( end != 0 ) end--; // we will now use 'end' as index into the array of stats

        if ( byLane[end] == null ) byLane[end] = new CycleStats(MAX_READ_LENGTH);
        if ( byLib[end] == null ) byLib[end] =new CycleStats(MAX_READ_LENGTH);
        byLane[end].add(quals);
        byLib[end].add(quals);

    }

    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {

        if ( AlignmentUtils.isReadUnmapped(read) && MAPPED_ONLY) return 0;

        SAMReadGroupRecord rg = read.getReadGroup();

        if ( rg == null ) throw new StingException("Read "+read.getReadName()+" is not assigned to any read group");

        String lane = read.getReadGroup().getPlatformUnit();
        String library = read.getReadGroup().getLibrary();

        if ( lane == null ) throw new StingException("Read "+read.getReadName()+" has no platform unit information");
        if ( library == null ) throw new StingException("Read "+read.getReadName()+" has no library information");

        int end = 0;

        if ( read.getReadPairedFlag() ) {

            if ( read.getFirstOfPairFlag() ) {
                if ( read.getSecondOfPairFlag() ) throw new StingException("Read "+read.getReadName()+" has conflicting first/second in pair attributes");
                end = 1;
            } else {
                if ( ! read.getSecondOfPairFlag() ) throw new StingException("Read "+read.getReadName()+" has conflicting first/second in pair attributes");
                end = 2;
            }
        }

        recordQuals(lane,library,AlignmentUtils.getQualsInCycleOrder(read),end,cyclesByLaneMap,cyclesByLibraryMap);

        return 1;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    public Integer reduceInit() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public Integer reduce(Integer value, Integer sum) {
        return sum.intValue()+value.intValue();  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void onTraversalDone(Integer result) {
        if ( HTML ) {
            out.println("<h3>Cycle Quality QC</h3>\n");
            out.println("File(s) analyzed: <br>");
            for ( File f : getToolkit().getArguments().samFiles) out.println(f.toString()+"<br>");
            out.println("<br>");
        }
        if ( HTML ) out.println("<br><br>");
        out.println("\n"+result+" reads analyzed\n");
        if ( HTML ) out.println("<br><br>");
        out.println("by platform unit:");
        if ( HTML ) out.println("<br>");
        report2(cyclesByLaneMap, new File(PREFIX+".byLane.txt"));
        out.println();
        if ( HTML ) out.println("<br><br>");
        out.println("\nby library:");
        if ( HTML ) out.println("<br>");
        report2(cyclesByLibraryMap, new File(PREFIX+".byLibrary.txt"));
        out.println();
        if ( HTML ) out.println("<br><br>");
    }

 

    private void report2(Map<String,CycleStats[]> m, File f) {
        long totalReads_1 =0;
        long totalReads_2 =0;
        long totalReads_unpaired = 0;
        SortedSet<String> columns = new TreeSet<String>();
        int maxLength = 0; // maximum read length across all lanes/read ends analyzed

        for( Map.Entry<String,CycleStats[]> e : m.entrySet() ) {
            if ( e.getValue()[0].getMaxReadLength() > maxLength ) maxLength = e.getValue()[0].getMaxReadLength();

            if ( e.getValue().length == 1 || e.getValue().length == 2 && e.getValue()[1] == null ) {
                totalReads_unpaired += e.getValue()[0].getReadCount(); // single end lane
            } else {
                totalReads_1 += e.getValue()[0].getReadCount();
                totalReads_2 += e.getValue()[1].getReadCount();
                if ( e.getValue()[1].getMaxReadLength() > maxLength ) maxLength = e.getValue()[1].getMaxReadLength();
            }

            columns.add(e.getKey());
        }

        if ( totalReads_1 == 0 && totalReads_2 != 0) {
            out.println("   End 1: No reads");
            if ( HTML ) out.println("<br>");
        }
        if ( totalReads_2 == 0 && totalReads_1 != 0 ) {
            out.println("   End 2: No reads");
            if ( HTML ) out.println("<br>");
        }
        if ( totalReads_1 == 0 && totalReads_2 == 0 && totalReads_unpaired == 0 ) {
            out.println("   No reads found.");
            if ( HTML ) out.println("<br>");
            return;
        }
        try {
            BufferedWriter w = new BufferedWriter(new FileWriter(f));

            w.write("cycle");

            for( String col : columns ) {
                CycleStats[] data = m.get(col);
                out.print("   ");
                out.print(col);

                CycleStats end1 = data[0];
                int minL = ( end1 == null ? 0 : end1.getMinReadLength() );
                int maxL = ( end1 == null ? 0 : end1.getMaxReadLength() );

                if ( data.length == 2 && data[1] != null ) {
                    out.println(": paired");
                    if ( HTML ) out.println("<br>");
                    out.println("    Reads analyzed:");
                    if ( HTML ) out.println("<br>");
                    CycleStats end2 = data[1];

                    out.print( "      End 1: "+ ( end1 == null ? 0 : end1.getReadCount()) );
                    if ( minL == maxL ) out.println("; read length = "+minL);
                    else out.println("; WARNING: variable read length = "+minL+"-"+maxL);
                    if ( HTML ) out.println("<br>");

                    out.print( "      End 2: "+ ( end2 == null ? 0 : end2.getReadCount()) );
                    minL = ( end2 == null ? 0 : end2.getMinReadLength() );
                    maxL = ( end2 == null ? 0 : end2.getMaxReadLength() );
                    if ( minL == maxL ) out.println("; read length = "+minL);
                    else out.println("; WARNING: variable read length = "+minL+"-"+maxL);
                    if ( HTML ) out.println("<br>");
                }
                else {
                    out.println(": unpaired");
                    if ( HTML ) out.println("<br>");
                    out.print( "      Reads analyzed: "+ ( end1 == null ? 0 : end1.getReadCount()) );
                    if ( minL == maxL ) out.println("; read length = "+minL);
                    else out.println("; WARNING: variable read length = "+minL+"-"+maxL);
                    if ( HTML ) out.println("<br>");
                }

                w.write('\t') ;
                w.write(col);
                if ( data.length == 1 ) {
                    w.write(".unpaired");
                }  else {
                    w.write(".end1");
                    w.write('\t') ;
                    w.write(col);
                    w.write(".end2");
                }
            }

            w.write('\n');

            int cycle = 0;

            Map<String,List<PrimitivePair.Int>> problems = new HashMap<String,List<PrimitivePair.Int>>();

            while ( cycle < maxLength ) {
                w.write(Integer.toString(cycle+1));
                for ( String col : columns ) {

                    CycleStats[] data = m.get(col);
                    CycleStats end1 = data[0];
                    w.write('\t');
                    if ( end1 == null || cycle >= end1.getMaxReadLength() ) w.write('.');
                    else {
                        double aq = end1.getAverageCycleQual(cycle);
                        w.write(String.format("%.4f",aq));
                        recordProblem(aq,cycle, problems,col+".End1");
                    }
                    if ( data.length > 1 && data[1] != null ) {
                        w.write('\t');
                        CycleStats end2 = data[1];
                        if ( end2 == null || cycle >= end2.getMaxReadLength() ) w.write('.');
                        else {
                            double aq = end2.getAverageCycleQual(cycle);
                            w.write(String.format("%.4f",aq));
                            recordProblem(aq,cycle, problems,col+".End2");
                        }
                    }
                }
                w.write('\n');
                cycle++;
            }
            w.close();

            if ( HTML ) out.println("<hr>");

            if ( HTML ) out.println("<br>");
            out.println("\nOUTCOME (threshold at Q="+QTHRESHOLD+"):");
            if ( HTML ) out.println("<br>");
            for ( String col : columns ) {
                List<PrimitivePair.Int> lp = problems.get(col+".End1");
                out.print("  "+col+" End1:");
                if ( lp == null ) {
                    out.print(" GOOD");
                } else {
                    for ( PrimitivePair.Int p : lp ) {
                        out.print(" "+(p.first+1)+"-");
                        if ( p.second >= 0 ) out.print((p.second+1));
                        else out.print("END");
                    }
                }
                out.println();
                if ( HTML ) out.println("<br>");

                lp = problems.get(col+".End2");
                out.print("  "+col+" End2:");
                if ( lp == null ) {
                    out.print(" GOOD");
                } else {
                    for ( PrimitivePair.Int p : lp ) {
                        out.print(" "+(p.first+1)+"-");
                        if ( p.second >= 0 ) out.print(p.second);
                        else out.print("END");
                    }
                }
                out.println();
                if ( HTML ) out.println("<br>");
            }

        } catch (IOException ioe) {
            throw new StingException("Failed to write report into the file "+f+":\n"+ioe.getMessage());
        }
    }


    private void recordProblem(double q, int cycle, Map<String,List<PrimitivePair.Int>> problems, String name) {

        PrimitivePair.Int p = null;
        List<PrimitivePair.Int> lp = null;
        if ( q < QTHRESHOLD ) { // there is a problem
               if ( ! problems.containsKey(name) ) {
                   lp = new ArrayList<PrimitivePair.Int>();
                   p = new PrimitivePair.Int(cycle,-1);
                   lp.add(p);
                   problems.put(name,lp);
               } else {
                   lp = problems.get(name);
                   p = lp.get(lp.size()-1);
               }
               if ( p.second != -1 ) { // if we are not already inside a run of bad qual bases
                   lp.add(new PrimitivePair.Int(cycle,-1)); // start new run
               }
        } else { // good base
              if ( problems.containsKey(name) ) { // only if we had problem intervals at all, we need to check if the last one needs to be closed
                  lp = problems.get(name);
                  p = lp.get(lp.size()-1);
                  if ( p.second == -1 ) p.second = cycle - 1;
              }
        }
    }


    static class CycleStats {
        private long readCount = 0;
        private long[] cycleQuals = null;
        private int minL = 1000000000; // read min. length
        private int maxL = 0; // read max. length

        public CycleStats(int N) {
            readCount = 0;
            cycleQuals = new long[N];
        }

        public void add(byte[] quals) {
            if ( quals.length > cycleQuals.length ) throw new StingException("A read of length "+quals.length+
                    " encountered, which exceeds specified maximum read length");
            if ( quals.length > maxL ) maxL = quals.length;
            if ( quals.length < minL ) minL = quals.length;
            for ( int i = 0 ; i < quals.length ; i++ ) {
                cycleQuals[i] += quals[i];
            }
            readCount++;
        }

        public long getReadCount() { return readCount; }
        public int getMaxReadLength() { return maxL; }
        public int getMinReadLength() { return minL; }
        long [] getCycleQualSums() { return cycleQuals; }
        long getCycleQualSum(int i) { return cycleQuals[i]; }
        double getAverageCycleQual(int i) { return ((double)cycleQuals[i])/readCount; }
    }
}
