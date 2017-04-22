/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.engine.recalibration.covariates;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.gatk.engine.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.gatk.engine.recalibration.ReadCovariates;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.collections.Pair;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public abstract class RepeatCovariate implements ExperimentalCovariate {
    protected int MAX_REPEAT_LENGTH;
    protected int MAX_STR_UNIT_LENGTH;
    private final HashMap<String, Integer> repeatLookupTable = new HashMap<String, Integer>();
    private final HashMap<Integer, String> repeatReverseLookupTable = new HashMap<Integer, String>();
    private int nextId = 0;

    // Initialize any member variables using the command-line arguments passed to the walkers
    @Override
    public void initialize(final RecalibrationArgumentCollection RAC) {
        MAX_STR_UNIT_LENGTH = RAC.MAX_STR_UNIT_LENGTH;
        MAX_REPEAT_LENGTH = RAC.MAX_REPEAT_LENGTH;
    }

    public void initialize(final int MAX_STR_UNIT_LENGTH, final int MAX_REPEAT_LENGTH) {
        this.MAX_STR_UNIT_LENGTH = MAX_STR_UNIT_LENGTH;
        this.MAX_REPEAT_LENGTH = MAX_REPEAT_LENGTH;
    }

    @Override
    public void recordValues(final GATKSAMRecord read, final ReadCovariates values) {
        // store the original bases and then write Ns over low quality ones
        final byte[] originalBases = read.getReadBases().clone();

        final boolean negativeStrand = read.getReadNegativeStrandFlag();
        byte[] bases = read.getReadBases();
        if (negativeStrand)
            bases = BaseUtils.simpleReverseComplement(bases);

        // don't record reads with N's
        if (!BaseUtils.isAllRegularBases(bases))
            return;

        for (int i = 0; i < bases.length; i++) {
            final Pair<byte[], Integer> res = findTandemRepeatUnits(bases, i);
            // to merge repeat unit and repeat length to get covariate value:
            final String repeatID =  getCovariateValueFromUnitAndLength(res.first,  res.second);
            final int key = keyForRepeat(repeatID);

            final int readOffset = (negativeStrand ? bases.length - i - 1 : i);
            values.addCovariate(key, key, key, readOffset);
        }

        // put the original bases back in
        read.setReadBases(originalBases);

    }

    /**
     * Please use {@link org.broadinstitute.gatk.utils.variant.TandemRepeatFinder#findMostRelevantTandemRepeatUnitAt(int)}
     * @deprecated
     */
    @Deprecated
    public Pair<byte[], Integer> findTandemRepeatUnits(byte[] readBases, int offset) {
        int maxBW = 0;
        byte[] bestBWRepeatUnit = new byte[]{readBases[offset]};
        for (int str = 1; str <= MAX_STR_UNIT_LENGTH; str++) {
            // fix repeat unit length
            //edge case: if candidate tandem repeat unit falls beyond edge of read, skip
            if (offset+1-str < 0)
                break;

            // get backward repeat unit and # repeats
            byte[] backwardRepeatUnit = Arrays.copyOfRange(readBases, offset - str + 1, offset + 1);
            maxBW = GATKVariantContextUtils.findNumberOfRepetitions(backwardRepeatUnit, Arrays.copyOfRange(readBases, 0, offset + 1), false);
            if (maxBW > 1) {
                bestBWRepeatUnit = backwardRepeatUnit.clone();
                break;
            }
        }
        byte[] bestRepeatUnit = bestBWRepeatUnit;
        int maxRL = maxBW;

        if (offset < readBases.length-1) {
            byte[] bestFWRepeatUnit = new byte[]{readBases[offset+1]};
            int maxFW = 0;
            for (int str = 1; str <= MAX_STR_UNIT_LENGTH; str++) {
                // fix repeat unit length
                //edge case: if candidate tandem repeat unit falls beyond edge of read, skip
                if (offset+str+1 > readBases.length)
                    break;

                // get forward repeat unit and # repeats
                byte[] forwardRepeatUnit = Arrays.copyOfRange(readBases, offset +1, offset+str+1);
                maxFW = GATKVariantContextUtils.findNumberOfRepetitions(forwardRepeatUnit, Arrays.copyOfRange(readBases, offset + 1, readBases.length), true);
                if (maxFW > 1) {
                    bestFWRepeatUnit = forwardRepeatUnit.clone();
                    break;
                }
            }
            // if FW repeat unit = BW repeat unit it means we're in the middle of a tandem repeat - add FW and BW components
            if (Arrays.equals(bestFWRepeatUnit, bestBWRepeatUnit)) {
                maxRL = maxBW + maxFW;
                bestRepeatUnit = bestFWRepeatUnit; // arbitrary
            }
            else {
                // tandem repeat starting forward from current offset.
                // It could be the case that best BW unit was differnet from FW unit, but that BW still contains FW unit.
                // For example, TTCTT(C) CCC - at (C) place, best BW unit is (TTC)2, best FW unit is (C)3.
                // but correct representation at that place might be (C)4.
                // Hence, if the FW and BW units don't match, check if BW unit can still be a part of FW unit and add
                // representations to total
                maxBW = GATKVariantContextUtils.findNumberOfRepetitions(bestFWRepeatUnit, Arrays.copyOfRange(readBases, 0, offset + 1), false);
                maxRL = maxFW + maxBW;
                bestRepeatUnit = bestFWRepeatUnit;

            }

        }



        if(maxRL > MAX_REPEAT_LENGTH) { maxRL = MAX_REPEAT_LENGTH; }
        return new Pair<byte[], Integer>(bestRepeatUnit, maxRL);

    }
    @Override
    public final Object getValue(final String str) {
        return str;
    }

    @Override
    public synchronized String formatKey(final int key) {
        // This method is synchronized so that we don't attempt to do a get()
        // from the reverse lookup table while that table is being updated
        return repeatReverseLookupTable.get(key);
    }

    @Requires({"repeatLength>=0", "repeatFromUnitAndLength != null"})
    @Ensures("result != null")
    protected abstract String getCovariateValueFromUnitAndLength(final byte[] repeatFromUnitAndLength, final int repeatLength);


    @Override
    public int keyFromValue(final Object value) {
        return keyForRepeat((String) value);
    }

    /**
     * Get the mapping from read group names to integer key values for all read groups in this covariate
     * @return a set of mappings from read group names -> integer key values
     */
    public Set<Map.Entry<String, Integer>> getKeyMap() {
        return repeatLookupTable.entrySet();
    }

    private int keyForRepeat(final String repeatID) {
        // Rather than synchronize this entire method (which would be VERY expensive for walkers like the BQSR),
        // synchronize only the table updates.

        // Before entering the synchronized block, check to see if this read group is not in our tables.
        // If it's not, either we will have to insert it, OR another thread will insert it first.
        // This preliminary check avoids doing any synchronization most of the time.
        if ( ! repeatLookupTable.containsKey(repeatID) ) {

            synchronized ( this ) {

                // Now we need to make sure the key is STILL not there, since another thread may have come along
                // and inserted it while we were waiting to enter this synchronized block!
                if ( ! repeatLookupTable.containsKey(repeatID) ) {
                    repeatLookupTable.put(repeatID, nextId);
                    repeatReverseLookupTable.put(nextId, repeatID);
                    nextId++;
                }
            }
        }

        return repeatLookupTable.get(repeatID);
    }


    /**
     * Splits repeat unit and num repetitions from covariate value.
     * For example, if value if "ATG4" it returns (ATG,4)
     * @param value             Covariate value
     * @return                  Split pair
     */
    @Requires("value != null")
    @Ensures({"result.first != null","result.second>=0"})
    public static Pair<String,Integer> getRUandNRfromCovariate(final String value) {

        int k = 0;
        for ( k=0; k < value.length(); k++ ) {
            if (!BaseUtils.isRegularBase(value.getBytes()[k]))
                break;
        }
        Integer nr = Integer.valueOf(value.substring(k,value.length())); // will throw NumberFormatException if format illegal
        if (k == value.length() || nr <= 0)
            throw new IllegalStateException("Covariate is not of form (Repeat Unit) + Integer");

        return new Pair<String,Integer>(value.substring(0,k), nr);
    }

    /**
     * Gets bases from tandem repeat representation (Repeat Unit),(Number of Repeats).
     * For example, (AGC),3 returns AGCAGCAGC
     * @param repeatUnit    Tandem repeat unit
     * @param numRepeats    Number of repeats
     * @return              Expanded String
     */
    @Requires({"numRepeats > 0","repeatUnit != null"})
    @Ensures("result != null")
    public static String getBasesFromRUandNR(final String repeatUnit, final int numRepeats) {
        final StringBuilder sb = new StringBuilder();

        for (int i=0; i < numRepeats; i++)
            sb.append(repeatUnit);

        return sb.toString();
    }

    // version given covariate key
    public static String  getBasesFromRUandNR(final String covariateValue) {
        Pair<String,Integer> pair = getRUandNRfromCovariate(covariateValue);
        return getBasesFromRUandNR(pair.getFirst(), pair.getSecond());
    }

    @Override
    public abstract int maximumKeyValue();



}
