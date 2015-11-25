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

package org.broadinstitute.gatk.engine.phonehome;

import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.simpleframework.xml.Element;
import org.simpleframework.xml.ElementList;

import java.util.ArrayList;
import java.util.List;

/**
 * A helper class for formatting in XML the throwable chain starting at e.
 */
class GATKRunReportException {
    @Element(required = false, name = "message")
    String message = null;

    @ElementList(required = false, name = "stacktrace")
    final List<String> stackTrace = new ArrayList<String>();

    @Element(required = false, name = "cause")
    GATKRunReportException cause = null;

    @Element(required = false, name = "is-user-exception")
    Boolean isUserException;

    @Element(required = false, name = "exception-class")
    Class exceptionClass;

    /**
     * Allow us to deserialize from XML
     */
    public GATKRunReportException() { }

    public GATKRunReportException(Throwable e) {
        message = e.getMessage();
        exceptionClass = e.getClass();
        isUserException = e instanceof UserException;
        for (StackTraceElement element : e.getStackTrace()) {
            stackTrace.add(element.toString());
        }

        if ( e.getCause() != null ) {
            cause = new GATKRunReportException(e.getCause());
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        GATKRunReportException that = (GATKRunReportException) o;

        if (cause != null ? !cause.equals(that.cause) : that.cause != null) return false;
        if (exceptionClass != null ? !exceptionClass.equals(that.exceptionClass) : that.exceptionClass != null)
            return false;
        if (isUserException != null ? !isUserException.equals(that.isUserException) : that.isUserException != null)
            return false;
        if (message != null ? !message.equals(that.message) : that.message != null) return false;
        if (stackTrace != null ? !stackTrace.equals(that.stackTrace) : that.stackTrace != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = message != null ? message.hashCode() : 0;
        result = 31 * result + (stackTrace != null ? stackTrace.hashCode() : 0);
        result = 31 * result + (cause != null ? cause.hashCode() : 0);
        result = 31 * result + (isUserException != null ? isUserException.hashCode() : 0);
        result = 31 * result + (exceptionClass != null ? exceptionClass.hashCode() : 0);
        return result;
    }
}
