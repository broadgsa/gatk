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

package org.broadinstitute.sting.oneoffprojects.utils;

import org.broadinstitute.sting.utils.exceptions.StingException;

import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;

public class AlignmentList implements Iterable<AlignmentInfo> {
        private int best_mm = 1000000000;
        private int next_best_mm = 1000000000;
        private List<AlignmentInfo> als = null;
        private int next_best_count = 0;
        private int best_overlap = 0;

        private AlignmentStrategy strategy = null;

        public AlignmentList(AlignmentStrategy s) {
            this.strategy = s;
            best_mm = 1000000000;
            next_best_mm = 1000000000;
            best_overlap = 0;
            als = new ArrayList<AlignmentInfo>(1);
        }

        public boolean isAligned() {
            return best_mm < 1000000000;
        }

        public List<AlignmentInfo> getAlignments() { return als; }

        public int size() { return als.size(); }

        public Iterator<AlignmentInfo> iterator() { return als.iterator(); }

 //       public void tryAdd(int mm, int offset, boolean isRc, int overlap) {
 //           tryAdd(new AlignmentInfo(mm,offset,isRc,overlap));
 //       }

        public void tryAdd(AlignmentInfo ai) {
            AlignmentStrategy.Action a = strategy.action(ai,this) ;
            switch ( a ) {
                case DISCARD: break;
                case REPLACE_BEST:
                    next_best_mm = best_mm;
                    next_best_count = size();
                    als.clear();
                    als.add(ai);
                    best_mm = ai.getMismatchCount();
                    best_overlap = ai.getOverlap();
                    break;
                case ADD_BEST:
                    als.add(ai);
                    if ( ai.getMismatchCount() < best_mm ) best_mm = ai.getMismatchCount();
                    if ( ai.getOverlap() > best_overlap) best_overlap = ai.getOverlap();
                    break;
                case REPLACE_NEXTBEST:
                    next_best_mm = ai.getMismatchCount();
                    next_best_count = 1;
                    break;
                case ADD_NEXTBEST:
                    next_best_count++;
                    if ( ai.getMismatchCount() < next_best_mm ) next_best_mm = ai.getMismatchCount();
                    break;
                default: throw new StingException("Unrecognized action requested: "+a);
            }
        }

        public void tryAddAll(AlignmentList al) {
            for( AlignmentInfo ai : al) {
                tryAdd(ai);
            }
        }

        public int getBestMMCount() { return best_mm; }
        public int getBestOverlap() { return best_overlap; }
        public int getBestHitCount() { return als.size() ; }
        public int getNextBestHitCount() { return next_best_count; }
        public int getNextBestMMCount() { return next_best_mm; }
//        public int getOverlap() { return overlap; }
//        public int getOffset() { return offset; }
//        public boolean isNegativeStrand() { return rc; }

//        public double getMismatchRate() { return isAligned() ? ((double)best_mm)/overlap : 1.0 ; }

}