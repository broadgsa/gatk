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

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.exceptions.StingException;

import java.util.List;
import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Sep 13, 2010
 * Time: 5:53:33 PM
 * To change this template use File | Settings | File Templates.
 */
public class AssemblyGraph {

    private List<Assembly> sources;

    public AssemblyGraph(Assembly a) {
        sources = new LinkedList<Assembly>();
        sources.add(a);
    }

    /** Initializes assembly from the single specified read, and sets this assembly as the root of this
     * assembly graph
     * @param r read; must be aligned, otherwise exception will be thrown
     * @param K index (Kmer) size of the assembly that will be initialized with the read r
     */
    public AssemblyGraph(SAMRecord r, int K) {
        if (AlignmentUtils.isReadUnmapped(r))
            throw new StingException("Can not initialize assembly graph with unaligned read");
        sources = new LinkedList<Assembly>();
        sources.add( new Assembly(K,r.getReadBases(),r.getReadName(), r.getAlignmentStart()) );
    }

    public void add(SAMRecord r) {
        
    }

}
