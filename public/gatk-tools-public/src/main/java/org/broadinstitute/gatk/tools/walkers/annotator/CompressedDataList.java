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

package org.broadinstitute.gatk.tools.walkers.annotator;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * A class to represent data as a list of <value,count> pairs.  For example, the list 2,2,2,2,2,2,3,4,4,4,5,5
 * would be compressed as 2,6,3,1,4,3,5,2. The compressed list should be sorted in ascending order by value.
 *
 * Created by gauthier on 9/25/15.
 */
public class CompressedDataList<T>  implements Iterable<T> {
    protected Map<T,Integer> valueCounts = new HashMap<>();

    public Map<T,Integer> getValueCounts(){
        return valueCounts;
    }

    public boolean isEmpty(){
        return valueCounts.isEmpty();
    }

    @Override
    public Iterator<T> iterator(){
        Iterator<T> it = new Iterator<T>() {
            private Iterator<T> keySetIterator = valueCounts.keySet().iterator();
            private T currentKey = valueCounts.isEmpty() ? null : keySetIterator.next();
            private int currentValueIndex = 0;
            private int currentValueSize = valueCounts.isEmpty() ? 0 : valueCounts.get(currentKey);

            @Override
            public boolean hasNext() {
                return !valueCounts.isEmpty() && (keySetIterator.hasNext() || currentValueIndex < currentValueSize);
            }

            @Override
            public T next() {
                T retKey = currentKey;
                currentValueIndex++;
                if(currentValueIndex==currentValueSize){
                    if(keySetIterator.hasNext()) {
                        currentKey = keySetIterator.next();
                        currentValueIndex = 0;
                        currentValueSize = valueCounts.get(currentKey);
                    }
                }
                return retKey;
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
        return it;
    }

    @Override
    public String toString(){
        String str = "";
        Object[] keys = valueCounts.keySet().toArray();
        Arrays.sort(keys);
        for (Object i: keys){
            if(!str.isEmpty())
                str+=",";
            str+=(i+","+valueCounts.get(i));
        }
        return str;
    }

    public void add(final T val){
        add(val, 1);
    }

    public void add(final T val, final int count){
        if(valueCounts.containsKey(val)){
            valueCounts.put(val, valueCounts.get(val)+count);
        }
        else
            valueCounts.put(val, count);

    }

    public void add(final CompressedDataList<T> obj){
        for(Map.Entry<T, Integer> pair : obj.getValueCounts().entrySet()){
            this.add(pair.getKey(),pair.getValue());
        }
    }

}
