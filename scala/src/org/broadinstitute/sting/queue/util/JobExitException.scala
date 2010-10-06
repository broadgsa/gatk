package org.broadinstitute.sting.queue.util

import org.broadinstitute.sting.queue.QException

/**
 * Captures the exit code and error text from a failed process.
 */
class JobExitException(var exitText: String, var commandLine: Array[String], var exitCode: Int, var stdErr: String)
        extends QException("%s%nCommand line:%n%s%nExit code: %s%nStandard error contained: %n%s"
                .format(exitText, commandLine.mkString(" "), exitCode, stdErr)) {
}
