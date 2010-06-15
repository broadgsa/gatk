package org.broadinstitute.sting.queue.util

import org.broadinstitute.sting.utils.text.XReadLines
import collection.mutable.ListBuffer
import collection.JavaConversions._
import java.io.File

object ProcessUtils extends Logging {

  Runtime.getRuntime.addShutdownHook(new Thread {
    override def run = for (process <- running.clone) {
      logger.warn("Killing: " + process)
      process.destroy
    }
  })

  val running = new ListBuffer[Process]() 

  def runCommandAndWait(command: String, directory: File) = {
    logger.debug("Running command: " + command)

    var builder = new ProcessBuilder("sh", "-c", command).directory(directory)

    var process = builder.start
    running += process
    var result = process.waitFor
    running -= process

    if (logger.isDebugEnabled) {
      for (line <- new XReadLines(process.getInputStream).iterator) {
        logger.debug("command: " + line)
      }

      for (line <- new XReadLines(process.getErrorStream).iterator) {
        logger.error("command: " + line)
      }
    }

    logger.debug("Command exited with result: " + result)

    result
  }
}
