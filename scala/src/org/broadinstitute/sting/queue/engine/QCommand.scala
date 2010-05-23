package org.broadinstitute.sting.queue.engine

/**
 * Defines a basic command to run
 * TODO: Allow overriding arguments per command such as the job queue
 */
class QCommand(val commandString: String) {
  override def toString = commandString
}
