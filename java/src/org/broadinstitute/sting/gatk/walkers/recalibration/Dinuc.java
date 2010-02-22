package org.broadinstitute.sting.gatk.walkers.recalibration;

/*
 * Copyright (c) 2009 The Broad Institute
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

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Nov 16, 2009
 */
public class Dinuc implements Comparable<Dinuc>{
    private byte first;
    private byte second;

    public Dinuc() {
        first = 0;
        second = 0;
    }

    public Dinuc(final byte _first, final byte _second) {
        first = _first;
        second = _second;
    }

    public final void setValues(final byte _first, final byte _second) {
        first = _first;
        second = _second;
    }

    public int compareTo(final Dinuc that) {
        if( this.first > that.first ) { return 1; }
        else if( this.first < that.first ) { return -1; }
        else { //this.first equals that.first
            if( this.second > that.second ) { return 1; }
            else if( this.second < that.second ) { return -1; }
            else { return 0; }
        }

    }

    public static int hashBytes(final byte byte1, final byte byte2) {
        return byte1 << 8 + byte2;
    }

    public String toString() { // This method call is how the Dinuc will get written out to the table recalibration file
        byte[] byteArray = {first,second};
        return new String(byteArray);
    }
}
