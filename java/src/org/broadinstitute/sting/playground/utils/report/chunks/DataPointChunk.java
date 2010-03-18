/*
 * Copyright (c) 2010.  The Broad Institute
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
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.utils.report.chunks;

import org.broadinstitute.sting.playground.utils.report.tags.DataPoint;
import org.broadinstitute.sting.utils.StingException;

import java.lang.reflect.Field;
import java.util.Map;


/**
 * @author aaron
 *         <p/>
 *         Class DataPointChunk
 *         <p/>
 *         a data point chunk
 */
public class DataPointChunk implements Chunk {
    public String name;
    public String description;
    public String value;

    /**
     * create a DataPointChunk from extracted information
     *
     * @param dpf           the datapoint field
     * @param mappings      the mapping of datapoints to fields in the analysis object
     * @param toMarshall    the object we're marshaling from
     */
    public DataPointChunk(Field dpf, Map<Field, DataPoint> mappings, Object toMarshall) {
        try {
            // make sure we can access the field
            dpf.setAccessible(true);
            name = mappings.get(dpf).name().equals("") ? dpf.getName() : mappings.get(dpf).name();
            description = mappings.get(dpf).description();
            value = dpf.get(toMarshall).toString();
        } catch (IllegalAccessException e) {
            throw new StingException("Unable to access variable " + mappings.get(mappings.get(dpf)), e);
        }
    }

    public String getName() {
        return name;
    }

    public String getDescription() {
        return description;
    }

    public String getValue() {
        return value;
    }

    /**
     * is the chunk we've created valid?  Invalid chunk would contain null data, or various other
     * factors that eliminate a chunk from being outputted.
     *
     * @return true if it's valid, false if not.
     */
    @Override
    public boolean isValid() {
        if (name == null || name.equals("")) return false;
        if (description == null || description.equals("")) return false;
        return value != null;
    }
}
