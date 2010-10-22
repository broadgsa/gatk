package org.broadinstitute.sting.queue.util

import org.broadinstitute.sting.queue.QException

/**
 * Captures the exit code and error text from a failed process.
 */
class JobExitException(val exitText: String, val commandLine: Array[String], val exitCode: Int, val stdErr: String)
        extends QException("%s%nCommand line:%n%s%nExit code: %s%nStandard error contained: %n%s"
                .format(exitText, commandLine.mkString(" "), exitCode, stdErr)) {
}
