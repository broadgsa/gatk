package org.broadinstitute.sting.queue.util

import org.broadinstitute.sting.queue.util.TextFormatUtils._

/**
 * Utilities for retrying a function and waiting a number of minutes.
 */
object Retry extends Logging {
  /**
   * Attempt the function and return the value on success.
   * Each time the function fails the stack trace is logged.
   * @param f Function to run and return its value.
   * @param wait The length in minutes to wait before each attempt.
   * @throws RetryException when all retries are exhausted.
   * @return The successful result of f.
   */
  def attempt[A](f: () => A, wait: Double*): A = {
    var count = 0
    var success = false
    val tries = wait.size + 1
    var result = null.asInstanceOf[A]
    while (!success && count < tries) {
      try {
        result = f()
        success = true
      } catch {
        case e => {
          count += 1
          if (count < tries) {
            val minutes = wait(count-1)
            logger.error("Caught error during attempt %s of %s."
                    .format(count, tries), e)
            logger.error("Retrying in %.1f minute%s.".format(minutes, plural(minutes)))
            Thread.sleep((minutes * 1000 * 60).toLong)
          } else {
            logger.error("Caught error during attempt %s of %s. Giving up."
                    .format(count, tries), e)
            throw new RetryException("Gave up after %s attempts.".format(tries), e)
          }
        }
      }
    }
    result
  }
}
