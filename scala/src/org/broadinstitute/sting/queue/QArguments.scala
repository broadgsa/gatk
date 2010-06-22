package org.broadinstitute.sting.queue

import collection.mutable.ListBuffer
import collection.JavaConversions._
import org.broadinstitute.sting.queue.util.Logging
import org.broadinstitute.sting.utils.text.XReadLines
import java.io.{FileInputStream, File}
import java.util.Properties

class QArguments(args: Array[String]) {
  var bsubAllJobs = false
  var dryRun = false
  val scripts = new ListBuffer[String]
  var inputPaths = List.empty[File]
  var argMap = Map.empty[String, String]

  val userArgs = parseArgs(args)

  private def parseArgs(args: Array[String]) = {
    var filtered = new ListBuffer[String]
    filtered.appendAll(args)

    if (isFlagged(filtered, "-debug"))
      Logging.setDebug
    if (isFlagged(filtered, "-trace"))
      Logging.setTrace
    if (isFlagged(filtered, "-dry"))
      dryRun = true
    if (isFlagged(filtered, "-bsub"))
      bsubAllJobs = true
    for (arg <- getArgs(filtered, "-P"))
      addArg(arg)
    for (arg <- getArgs(filtered, "-I"))
      addFile(arg)
    for (arg <- getArgs(filtered, "-S"))
      scripts.append(arg)

    List(filtered:_*)
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
      argMap += tokens(0) -> tokens(1)
    } else if (arg.endsWith(".properties")) {
      if (!file.exists)
        throw new QException("File not found: " + file.getAbsolutePath)
      var props = new Properties
      props.load(new FileInputStream(file))
      for ((name, value) <- props)
        argMap += name -> value
    } else {
      throw new QException("Invalid property: " + arg)
    }
  }

  def addFile(arg: String): Unit = {
    var file = new File(arg)
    inputPaths :+= file
    if (arg.endsWith(".list"))
      new XReadLines(file).iterator.foreach(addFile(_))
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
