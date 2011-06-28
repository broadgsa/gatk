package org.broadinstitute.sting.queue.util

object TextFormatUtils {
  /**
   * Returns the string "s" if x is greater than 1.
   * @param x Value to test.
   * @return "s" if x is greater than one else "".
   */
  def plural(x: Int) = if (x > 1) "s" else ""

  /**
   * Returns the string "s" if x is greater than 1.
   * @param x Value to test.
   * @return "s" if x is greater than one else "".
   */
  def plural(x: Long) = if (x > 1) "s" else ""

  /**
   * Returns the string "s" if x is not equal to 1.
   * @param x Value to test.
   * @return "s" if x is greater than one else "".
   */
  def plural(x: Float) = if (x != 1) "s" else ""

  /**
   * Returns the string "s" if x is not equal to 1.
   * @param x Value to test.
   * @return "s" if x is greater than one else "".
   */
  def plural(x: Double) = if (x != 1) "s" else ""
}
