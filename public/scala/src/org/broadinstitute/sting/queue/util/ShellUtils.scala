package org.broadinstitute.sting.queue.util

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