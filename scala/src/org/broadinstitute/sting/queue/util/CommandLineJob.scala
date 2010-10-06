package org.broadinstitute.sting.queue.util

import java.io.File

/**
 * Base utility class for a command line job.
 */
abstract class CommandLineJob {
  var command: String = _
  var workingDir: File = _
  var inputFile: File = _
  var outputFile: File = _
  var errorFile: File = _

  /**
   * Runs the command, either immediately or dispatching it to a compute farm.
   * If it is dispatched to a compute farm it should not start until jobs it depends on are finished.
   */
  def run()

  /**
   * Returns the content of a command output.
   * @param streamOutput The output of the command.
   * @return The content of the command, along with a message if it was truncated.
   */
  protected def content(streamOutput: ProcessController.StreamOutput) = {
    var content = streamOutput.content
    if (streamOutput.contentTruncated)
      content += "%n%n<truncated>".format()
    content
  }

  /**
   * Returns the ProcessController for this thread.
   * @return The ProcessController for this thread.
   */
  protected def processController = CommandLineJob.threadProcessController.get

  /** A five mb limit of characters for display. */
  protected val FIVE_MB = 1024 * 512 * 5;
}

/**
 * Base class for a command line job.
 */
object CommandLineJob {
  /** Thread local process controller container. */
  private val threadProcessController = new ThreadLocal[ProcessController] {
    override def initialValue = new ProcessController
  }
}
