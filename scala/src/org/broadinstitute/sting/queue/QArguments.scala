package org.broadinstitute.sting.queue

import collection.mutable.ListBuffer
import collection.JavaConversions._
import org.broadinstitute.sting.queue.util.Logging
import org.broadinstitute.sting.utils.text.XReadLines
import java.io.{FileInputStream, File}
import java.util.Properties

class QArguments(args: Array[String]) {
  var useBsub = false
  var dryRun = false
  val scripts = new ListBuffer[String]
  var inputPaths = List.empty[File]
  var argMap = Map.empty[String, String]

  parseArgs(args)

  private def parseArgs(args: Array[String]) = {
    var filtered = new ListBuffer[String]
    filtered.appendAll(args)

    if (isFlagged(filtered, "-debug"))
      Logging.enableDebug

    if (isFlagged(filtered, "-dry"))
      dryRun = true
    if (isFlagged(filtered, "-bsub"))
      useBsub = true
    for (arg <- getArgs(filtered, "-P"))
      addArg(arg)
    for (arg <- getArgs(filtered, "-I"))
      addFile(arg)
    for (arg <- getArgs(filtered, "-S"))
      scripts.append(arg)
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

  def addArg(arg: String) = {
    var file = new File(arg)
    if (arg.contains("=") && !file.exists) {
      val tokens = arg.split("=", 2)
      argMap = argMap.updated(tokens(0), tokens(1))
    } else if (file.exists && arg.endsWith(".properties")) {
      var props = new Properties
      props.load(new FileInputStream(file))
      for ((name, value) <- props)
        argMap = argMap.updated(name, value)
    }
  }

  def addFile(arg: String): Unit = {
    var file = new File(arg)
    if (arg.endsWith(".list")) {
      new XReadLines(file).iterator.foreach(addFile(_))
    } else {
      inputPaths = inputPaths ::: List(file)
    }
  }
}

object QArguments {
  def strip(filtered: ListBuffer[String], search: String) = {
    var index = 0
    while (0 <= index && index < filtered.size) {
      index = filtered.indexOf(search)
      if (index >= 0) {
        filtered.remove(index, 2)
      }
    }
  }
}
