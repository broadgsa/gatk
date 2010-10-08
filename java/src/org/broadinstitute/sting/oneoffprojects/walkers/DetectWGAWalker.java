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

package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceDataSource;
import org.broadinstitute.sting.oneoffprojects.utils.ReadPair;
import org.broadinstitute.sting.oneoffprojects.utils.AlignmentInfo;
import org.broadinstitute.sting.oneoffprojects.utils.Assembly;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import net.sf.samtools.SAMRecord;

import java.util.Map;
import java.util.HashMap;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Sep 13, 2010
 * Time: 12:45:42 PM
 * To change this template use File | Settings | File Templates.
 */

@WalkerName("DetectWGA")
@Requires(value={DataSource.REFERENCE_BASES})

public class DetectWGAWalker extends ReadWalker<Integer,Integer> {

    private int TIP_LENGTH = 10;
    private int TIP_MM_THRESHOLD = 1;
    private double TIP_AV_QUAL_THRESHOLD = 15.0;
    
    private boolean DEBUG = true;

    Map<String, ReadPair> pairCache = null;
    Map<String,MathUtils.RunningAverage> fragmentSizeMap = null;   // by library
    private ReferenceDataSource refData;
    private byte[] refBases;


   

    @Override
    public void initialize() {
        refData = new ReferenceDataSource(getToolkit().getArguments().referenceFile);
        pairCache = new HashMap<String, ReadPair>();
        fragmentSizeMap = new HashMap<String,MathUtils.RunningAverage>();
    }


    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {

        if ( ! read.getReadPairedFlag() ) return 0; // for now!!

        // read is paired

        cacheReadAsPair(read);

        /*
        // if the read is already mapped (uniquely), we check if it may have the "colored tip" artifact on either side:
        if ( AlignmentUtils.isReadUniquelyMapped(read) ) {

            TipInfo tips = countTipMismatches(read,TIP_LENGTH);

            if ( tips.leftMM() >= TIP_MM_THRESHOLD || tips.rightMM() >= TIP_MM_THRESHOLD ) {
                if ( DEBUG ) {
                    out.println("  Read "+read.getReadName()+ " has "+tips.leftMM()+"/"+tips.rightMM()+" mismatches in the tips");
                    out.println("  Pair orientation: "+pair.getPairType());
                }
                // try adding read to existing assemblies:
                AlignmentInfo al = alignToAllAddToBest(read,Math.min(3,tips.leftMM()+tips.rightMM())-1);
                if ( al == null ) {
                    if ( tips.leftMM() >= TIP_MM_THRESHOLD && tips.leftQ() >= TIP_AV_QUAL_THRESHOLD ||
                        tips.rightMM() >= TIP_MM_THRESHOLD && tips.rightQ() >= TIP_AV_QUAL_THRESHOLD ) {
                        if ( DEBUG ) out.println("   Initialized new assembly.") ;
                        Assembly a = new Assembly(read.getReadBases(),read.getReadName(),read.getAlignmentStart());
                        tryAndAddUnmapped(a); // see if we got unmapped reads that would align nicely
                        assemblies.add(a);
                    }
                }
            }
            return 1;
        }
*/



        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    public Integer reduceInit() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public Integer reduce(Integer value, Integer sum) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /** little helper: if we already cached the pair object for this read, just add the read to that object; if we did not - instantiate
     * new pair objetc first and register it in the map, then add the read; this method also updates other cache(s)/trackers as needed, e.g.
     * fragment size map
     */
    private void cacheReadAsPair(SAMRecord read) {
        ReadPair pair = pairCache.get( read.getReadName() );
        if ( pair == null ) {
            pair = new ReadPair(read);
            pairCache.put(read.getReadName(),pair);
        }

        pair.addRead(read);

        // if it's a good pair, add its fragment size to the stats:
        if ( pair.hasBothEnds() && pair.bothEndsMapped() && pair.isProper() ) {
            String lib = read.getReadGroup().getLibrary();
            MathUtils.RunningAverage fSize = fragmentSizeMap.get(lib);
            if ( fSize == null ) {
                fSize = new MathUtils.RunningAverage();
                fragmentSizeMap.put(lib,fSize);
            }
            fSize.add(pair.getFragmentSize());
        }

    }

    private TipInfo countTipMismatches(SAMRecord read, int tip_length) {

            AlignmentUtils.MismatchCount left_mm = AlignmentUtils.getMismatchCount(read,refBases,read.getAlignmentStart()-1,0,tip_length);

            int right_start = read.getReadLength()-tip_length;
            AlignmentUtils.MismatchCount right_mm = AlignmentUtils.getMismatchCount(read,refBases,read.getAlignmentStart()-1,right_start,read.getReadLength()-right_start);

            return new TipInfo(left_mm,right_mm);
    }

    class TipInfo {
        AlignmentUtils.MismatchCount left_mm;
        AlignmentUtils.MismatchCount right_mm;
        double left_avQ;
        double right_avQ;

        public TipInfo(AlignmentUtils.MismatchCount l,AlignmentUtils.MismatchCount r) {
            left_mm = l;
            right_mm = r;
            left_avQ = (l.numMismatches ==0 ? 0 : ((double)l.mismatchQualities)/l.numMismatches );
            right_avQ = (r.numMismatches ==0 ? 0 : ((double)r.mismatchQualities)/r.numMismatches );
        }

        public int leftMM() { return left_mm.numMismatches; }
        public int rightMM() { return right_mm.numMismatches; }
        public double leftQ() { return left_avQ; }
        public double rightQ() { return right_avQ; }
    }

}
