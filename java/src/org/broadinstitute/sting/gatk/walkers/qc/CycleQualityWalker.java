/*
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the ”Software”), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED ”AS IS”, WITHOUT WARRANTY OF ANY KIND,
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
import org.broadinstitute.sting.utils.StingException;
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
        if ( byLib == null ) libMap.put(library, byLib =new CycleStats[(end==0?1:2)]);

        if ( end != 0 ) end--; // we will now use 'end' as index into the array of stats

        if ( byLane[end] == null ) byLane[end] = new CycleStats(MAX_READ_LENGTH);
        if ( byLib[end] == null ) byLib[end] =new CycleStats(MAX_READ_LENGTH);
        byLane[end].add(quals);
        byLib[end].add(quals);

    }

    public Integer map(char[] ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {

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
        out.println(result+" reads analyzed");
        out.println("by platform unit:");
        report2(cyclesByLaneMap, new File(PREFIX+".byLane.txt"));
        out.println("\nby library:");
        report2(cyclesByLibraryMap, new File(PREFIX+".byLibrary.txt"));
    }

    private void report(Map<String,CycleStats> m, File f) {
        long totalReads =0;
        for( Map.Entry<String,CycleStats> e : m.entrySet() ) {
            totalReads += e.getValue().getReadCount();
        }
        if ( totalReads == 0 ) {
            out.println("   No reads");
            return;
        }
        try {
            FileWriter w = new FileWriter(f);

            SortedSet<String> columns = new TreeSet<String>();

            int maxLength = 0;

            for( Map.Entry<String,CycleStats> e : m.entrySet() ) {
                out.print( "   "+ e.getKey()+ ": "+ e.getValue().getReadCount());
                int minL = e.getValue().getMinReadLength();
                int maxL = e.getValue().getMaxReadLength();
                if ( minL == maxL ) out.println("; read length = "+minL);
                else out.println("; WARNING: variable read length = "+minL+"-"+maxL);
                columns.add(e.getKey());
                if ( e.getValue().getMaxReadLength() > maxLength ) maxLength = e.getValue().getMaxReadLength();
            }

            w.write("cycle");
            for ( String col : columns ) { w.write('\t') ; w.write(col); }
            w.write('\n');

            int cycle = 0;

            while ( cycle < maxLength ) {
                w.write(Integer.toString(cycle+1));
                for ( String col : columns ) {
                    CycleStats cs = m.get(col);
                    w.write('\t');
                    if ( cycle >= cs.getMaxReadLength() ) w.write('.');
                    else w.write(String.format("%.4f",cs.getAverageCycleQual(cycle)));
                }
                w.write('\n');
                cycle++;
            }
            w.close();
        } catch (IOException ioe) {
            throw new StingException("Failed to write report into the file "+f+":\n"+ioe.getMessage());
        }
    }


    private void report2(Map<String,CycleStats[]> m, File f) {
        long totalReads_1 =0;
        long totalReads_2 =0;
        long totalReads_unpaired = 0;
        SortedSet<String> columns = new TreeSet<String>();
        int maxLength = 0; // maximum read length across all lanes/read ends analyzed

        for( Map.Entry<String,CycleStats[]> e : m.entrySet() ) {
            if ( e.getValue()[0].getMaxReadLength() > maxLength ) maxLength = e.getValue()[0].getMaxReadLength();

            if ( e.getValue().length == 1 ) {
                totalReads_unpaired += e.getValue()[0].getReadCount(); // single end lane
            } else {
                totalReads_1 += e.getValue()[0].getReadCount();
                totalReads_2 += e.getValue()[1].getReadCount();
                if ( e.getValue()[1].getMaxReadLength() > maxLength ) maxLength = e.getValue()[1].getMaxReadLength();
            }

            columns.add(e.getKey());
        }

        if ( totalReads_1 == 0 ) {
            out.println("   End 1: No reads");
        }
        if ( totalReads_2 == 0 ) {
            out.println("   End 2: No reads");
        }
        if ( totalReads_1 == 0 && totalReads_2 == 0 ) return;
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

                if ( data.length == 1 ) {
                    out.println(": paired");
                    CycleStats end2 = data[1];

                    out.print( "      End 1: "+ ( end1 == null ? 0 : end1.getReadCount()) );
                    if ( minL == maxL ) out.println("; read length = "+minL);
                    else out.println("; WARNING: variable read length = "+minL+"-"+maxL);

                    out.print( "      End 2: "+ ( end2 == null ? 0 : end2.getReadCount()) );
                    minL = ( end2 == null ? 0 : end2.getMinReadLength() );
                    maxL = ( end2 == null ? 0 : end2.getMaxReadLength() );
                    if ( minL == maxL ) out.println("; read length = "+minL);
                    else out.println("; WARNING: variable read length = "+minL+"-"+maxL);
                }
                else {
                    out.println(": unpaired");
                    out.print( "      Reads: "+ ( end1 == null ? 0 : end1.getReadCount()) );
                    if ( minL == maxL ) out.println("; read length = "+minL);
                    else out.println("; WARNING: variable read length = "+minL+"-"+maxL);
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

            while ( cycle < maxLength ) {
                w.write(Integer.toString(cycle+1));
                for ( String col : columns ) {

                    CycleStats[] data = m.get(col);
                    CycleStats end1 = data[0];
                    w.write('\t');
                    if ( end1 == null || cycle >= end1.getMaxReadLength() ) w.write('.');
                    else w.write(String.format("%.4f",end1.getAverageCycleQual(cycle)));
                    if ( data.length > 1 ) {
                        w.write('\t');
                        CycleStats end2 = data[1];
                        if ( end2 == null || cycle >= end2.getMaxReadLength() ) w.write('.');
                        else w.write(String.format("%.4f",end2.getAverageCycleQual(cycle)));
                    }
                }
                w.write('\n');
                cycle++;
            }
            w.close();
        } catch (IOException ioe) {
            throw new StingException("Failed to write report into the file "+f+":\n"+ioe.getMessage());
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
