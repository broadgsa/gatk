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

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Aug 3, 2010
 * Time: 4:53:01 PM
 * To change this template use File | Settings | File Templates.
 */
public class DefaultAlignmentStrategy implements AlignmentStrategy {

    public Action action(AlignmentInfo alignment, AlignmentList currentList) {
        if ( currentList.size() == 0 ) return Action.REPLACE_BEST;

        if ( alignment.getMismatchCount() > currentList.getNextBestMMCount() ) return Action.DISCARD;

        if ( alignment.getMismatchCount() < currentList.getBestMMCount() ) return Action.REPLACE_BEST;
        if ( alignment.getMismatchCount() == currentList.getBestMMCount() ) return Action.ADD_BEST;
        if ( alignment.getMismatchCount() < currentList.getNextBestMMCount() ) return Action.REPLACE_NEXTBEST;
        if ( alignment.getMismatchCount() == currentList.getNextBestMMCount() ) return Action.ADD_NEXTBEST;

        throw new StingException("Unexpected case found and left unprocessed");
//        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
