/*
 * Copyright (c) 2010 The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the ”Software”), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED ”AS IS”, WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.utils;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.util.CloseableIterator;

import java.util.List;
import java.util.ArrayList;

import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.gatk.iterators.PushbackIterator;

/**
 * Iterates synchronously over two SAM files. At each iteration returs alignments with the same read name (in the order
 * the read names appear in the files). Alignment(s) from the first/second SAM file will be stored as the first/second
 * element of the pair, respectively. Multiple alignments (alternative placements) for a given read are allowed in both
 * input files. If only one of the files have alignment(s) for a given read name, the returned
 * pair will contain an empty list in the element corresponding to the other file. To enable this sort of traversal
 * synchronized by read names, the input SAM files must be sorted by read name. Constructor of this class verifies
 * that this is the case: SAM file  headers must report either 'queryname' sorting order, or no sorting order
 * (to allow the code to work with not-fully compliant 3rd party tools that do not set header flags properly; a warning
 * will be currently printed to stdout in this case); if sorting order in either of the files is set to "coordinate",
 * an exception will be thrown.
 */
    public class ParallelSAMIterator implements CloseableIterator< Pair< List<SAMRecord>, List<SAMRecord> > > {
        private SAMFileReader reader1;
        private SAMFileReader reader2;
        PushbackIterator<SAMRecord> i1;
        PushbackIterator<SAMRecord> i2;
        List<SAMRecord> alignments1;
        List<SAMRecord> alignments2;

        public ParallelSAMIterator(SAMFileReader r1, SAMFileReader r2) {
            reader1 = r1;
            reader2 = r2;
            checkSortOrder(r1,"End 1");
            checkSortOrder(r2, "End 2");
            i1 = new PushbackIterator(r1.iterator());
            i2 = new PushbackIterator(r2.iterator());
            alignments1 = nextGroup(i1); // pre-read next set of alignments
            alignments2 = nextGroup(i2);
        }


        /**
         * Returns <tt>true</tt> if the iteration has more elements. (In other
         * words, returns <tt>true</tt> if <tt>next</tt> would return an element
         * rather than throwing an exception.)
         *
         * @return <tt>true</tt> if the iterator has more elements.
         */
        public boolean hasNext() {
            return alignments1.size() > 0  || alignments2.size() > 0;
        }

        /**
         * Returns the next element in the iteration.
         *
         * @return the next element in the iteration.
         * @throws java.util.NoSuchElementException
         *          iteration has no more elements.
         */
        public Pair< List<SAMRecord>, List<SAMRecord> > next() {
            Pair< List<SAMRecord>, List<SAMRecord> > result;

            if ( alignments1.size() == 0 ) {
                // no more alignments left for end1
                result =  new Pair< List<SAMRecord>, List<SAMRecord> >(alignments1,alignments2);
                alignments2 = nextGroup(i2);
                return result;
            }
            if ( alignments2.size() == 0 ) {
                // no more alignments left for end2
                result =  new Pair< List<SAMRecord>, List<SAMRecord> >(alignments1,alignments2);
                alignments1 = nextGroup(i1);
                return result;
            }
            // next group of alignments is held for both ends. Check the read names:
            String end1Name = alignments1.get(0).getReadName();
            String end2Name = alignments2.get(0).getReadName();

            int cmp = end1Name.compareTo(end2Name);
            if ( cmp < 0 ) {
                // end1 goes before end2; return end1 with empty list for corresponding end2 and read next end1
                result = new Pair< List<SAMRecord>, List<SAMRecord> >(alignments1,new ArrayList<SAMRecord>());
                alignments1 = nextGroup(i1);
            } else {
                if ( cmp > 0 ) {
                    // end2 goes before end1; return end2 with empty list for corresponding end1 and read next end2
                    result = new Pair< List<SAMRecord>, List<SAMRecord> >(new ArrayList<SAMRecord>(),alignments2);
                    alignments2 = nextGroup(i2);
                } else {
                    // end 1 and end2 have the same read name => we got a mate pair:
                    result = new Pair< List<SAMRecord>, List<SAMRecord> >(alignments1, alignments2);
                    alignments1 = nextGroup(i1);
                    alignments2 = nextGroup(i2);
                }
            }
            return result;
        }

        /**
         * Removes from the underlying collection the last element returned by the
         * iterator (optional operation).  This method can be called only once per
         * call to <tt>next</tt>.  The behavior of an iterator is unspecified if
         * the underlying collection is modified while the iteration is in
         * progress in any way other than by calling this method.
         *
         * @throws UnsupportedOperationException if the <tt>remove</tt>
         *                                       operation is not supported by this Iterator.
         * @throws IllegalStateException         if the <tt>next</tt> method has not
         *                                       yet been called, or the <tt>remove</tt> method has already
         *                                       been called after the last call to the <tt>next</tt>
         *                                       method.
         */
        public void remove() {
            throw new UnsupportedOperationException("ParallelSAMIterator does not support remove() operation.");
        }

        public void close() {
            reader1.close();
            reader2.close();
        }

        /**
         * Read next alignment, and all immediately following ones that share same read name with the first;
         * return them all as a list.
         * @param i
         * @return
         */
        private List<SAMRecord> nextGroup(PushbackIterator<SAMRecord> i) {
            List<SAMRecord> result = new ArrayList<SAMRecord>();
            String readName ;

            if ( ! i.hasNext() ) return result; // nothing left
            SAMRecord r = i.next();
            readName = r.getReadName();
            result.add(r);

            while ( i.hasNext() ) {
                r = i.next();
                if ( ! r.getReadName().equals(readName) ) {
                    i.pushback(r);
                    break;
                }
                result.add(r);
            }
            return result;
        }

        /**
         * Utility method: checks that the sorting order in the input file is right
        *
        * @param reader sam file reader
        * @param fileName name of the file the reader is associated with. Used only to create more intelligible warning/exception messages,
        * you can actually pass any string here.
        */
        private void checkSortOrder(SAMFileReader reader, String fileName) {

            if  ( reader.getFileHeader() == null ) {
                System.out.println("WARNING: File "+fileName+" has no header. Assuming that file is sorted by read name.");
            }

            switch ( reader.getFileHeader().getSortOrder() ) {
            case coordinate:
               throw new RuntimeException("File "+fileName+" is sorted by coordinate. Sort it by read name first.");
            case unsorted:
               System.out.println("WARNING: file "+fileName+" has sorting order tag set to 'unsorted'. "+
                        "Assuming that it is sorted by read name.");
                break;
            case queryname: break; // good, that's what we need
                 default: throw new RuntimeException("File "+fileName + ": unknown sorting order ("+
                        reader.getFileHeader().getSortOrder()+")");
            }

        }

    }
