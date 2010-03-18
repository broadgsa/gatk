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

import org.broadinstitute.sting.playground.utils.report.TableFormatter;
import org.broadinstitute.sting.playground.utils.report.tags.Table;
import org.broadinstitute.sting.utils.StingException;

import java.lang.reflect.Field;
import java.util.List;
import java.util.Map;


/**
 * @author aaron
 *         <p/>
 *         Class TableChunk
 *         <p/>
 *         all the information we gleam from a Table annotation, including the table entries
 */
public class TableChunk implements Chunk {
    public final String name;
    public final String header;
    public final String description;
    public final int columns;
    public final TableFormatter table;

    /**
     * create a TableChunk from extracted information
     *
     * @param tbf         the Table annotated field
     * @param mappings   the mapping of datapoints to fields in the analysis object
     * @param toMarshall the object we're marshaling from
     */
    public TableChunk(Field tbf, Map<Field, Table> mappings, Object toMarshall) {
        try {
            // make sure we can access the field
            tbf.setAccessible(true);
            name = mappings.get(tbf).name().equals("") ? tbf.getName() : mappings.get(tbf).name();
            description = mappings.get(tbf).description();
            header = mappings.get(tbf).header();
            columns = mappings.get(tbf).columns();
            table = new TableFormatter(tbf.get(toMarshall));
        } catch (IllegalAccessException e) {
            throw new StingException("Unable to access variable " + tbf, e);
        }
    }

    public String getName() {
        return name;
    }

    public String getHeader() {
        return header;
    }

    public String getDescription() {
        return description;
    }

    public int getColumns() {
        return columns;
    }

    public TableFormatter getTable() {
        return table;
    }

    public List<List<Object>> getRows() {
        return table.toRows();    
    }

    /**
     * is the chunk we've created valid?  Invalid chunk would contain null data, or various other
     * factors that eliminate a chunk from being outputted.
     *
     * @return true if it's valid, false if not.
     */
    @Override
    public boolean isValid() {
        if (this.name == null || this.name.equals(""))      return false;
        if (this.header == null || this.header.equals(""))  return false;
        if (this.description == null)                       return false;
        if (columns < 1)                                    return false;
        if (table == null || table.tableEntries == null)    return false;
        return true;
    }
}
