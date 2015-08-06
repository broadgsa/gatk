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

package org.broadinstitute.gatk.utils.exceptions;

import java.lang.reflect.InvocationTargetException;

/**
 * Class for handling common failures of dynamic class resolution
 */
public class DynamicClassResolutionException extends UserException {
    public DynamicClassResolutionException(Class c, Exception ex) {
        super(String.format("Could not create module %s because %s caused by exception %s",
                c.getSimpleName(), moreInfo(ex), ex.getMessage()));
    }

    private static String moreInfo(Exception ex) {
        try {
            throw ex;
        } catch (InstantiationException e) {
            return "BUG: cannot instantiate class: must be concrete class";
        } catch (NoSuchMethodException e) {
            return "BUG: Cannot find expected constructor for class";
        } catch (IllegalAccessException e) {
            return "Cannot instantiate class (Illegal Access)";
        } catch (InvocationTargetException e) {
            return "Cannot instantiate class (Invocation failure)";
        } catch ( Exception e ) {
            return String.format("an exception of type %s occurred",e.getClass().getSimpleName());
        }
    }
}
