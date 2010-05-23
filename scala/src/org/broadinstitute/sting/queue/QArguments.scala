package org.broadinstitute.sting.queue

import collection.mutable.ListBuffer
import engine.Pipeline
import org.broadinstitute.sting.queue.util.Logging

object QArguments {
  val scripts = new ListBuffer[String]

  def parseArgs(args: Array[String]) {
    filterArgs(args)
  }

  /**
   * Pull out any args that are meant for QCommandLine
   */
  private def filterArgs(args: Array[String]) = {
    var filtered = new ListBuffer[String]
    filtered.appendAll(args)

    if (isFlagged(filtered, "-debug"))
      Logging.enableDebug

    if (isFlagged(filtered, "-dry"))
      Pipeline.dryRun = true
    if (isFlagged(filtered, "-bsub"))
      Pipeline.useBsub = true
    for (arg <- getArgs(filtered, "-P"))
      Pipeline.addArg(arg)
    for (arg <- getArgs(filtered, "-I"))
      Pipeline.addFile(arg)
    for (arg <- getArgs(filtered, "-S"))
      scripts.append(arg)

    filtered
  }

  private def isFlagged(filtered: ListBuffer[String], search: String) = {
    var found = false
    var index = 0
    while (0 <= index && index < filtered.size) {
      index = filtered.indexOf(search)
      if (index >= 0) {
        found = true
        filtered.remove(index)
      }
    }
    found
  }

  private def getArgs(filtered: ListBuffer[String], search: String) = {
    var found = new ListBuffer[String]
    var index = 0
    while (0 <= index && index < filtered.size) {
      index = filtered.indexOf(search)
      if (index >= 0) {
        found.append(filtered(index+1))
        filtered.remove(index, 2)
      }
    }
    found
  }


  def strip(filtered: ListBuffer[String], search: String) = {
    var index = 0
    while (0 <= index && index < filtered.size) {
      index = filtered.indexOf(search)
      if (index >= 0) {
        filtered.remove(index, 2)
      }
    }
  }}
