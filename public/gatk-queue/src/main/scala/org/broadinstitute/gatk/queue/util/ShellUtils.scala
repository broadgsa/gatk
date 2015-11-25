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

package org.broadinstitute.gatk.queue.util

import java.lang.IllegalArgumentException

object ShellUtils {

  /**
   * Escapes the String it's passed so that it will be interpreted literally when
   * parsed by sh/bash. Can correctly escape all characters except \0, \r, and \n
   *
   * Replaces all instances of ' with '\'', and then surrounds the resulting String
   * with single quotes.
   *
   * Examples:
   * ab   ->  'ab'
   * a'b  ->  'a'\''b'
   * ''   ->  ''\'''\'''
   *
   * Since \' is not supported inside single quotes in the shell (ie., '\'' does not work),
   * whenever we encounter a single quote we need to terminate the existing single-quoted
   * string, place the \' outside of single quotes, and then start a new single-quoted
   * string. As long as we don't insert spaces between the separate strings, the shell will
   * concatenate them together into a single argument value for us.
   *
   * @param str the String to escape
   * @return the same String quoted so that it will be interpreted literally when
   *         parsed by sh/bash
   */
  def escapeShellArgument ( str : String ) : String = {
    if ( str == null ) {
      throw new IllegalArgumentException("escapeShellArgument() was passed a null String")
    }

    "'" + str.replaceAll("'", "'\\\\''") + "'"
  }
}