package org.broadinstitute.sting.queue.engine

import graphing.JobGrapher
import org.broadinstitute.sting.queue.util.Logging
import scheduling.{DispatchJobScheduler, SimpleJobScheduler}
import java.lang.String
import collection.immutable.List
import collection.JavaConversions._
import java.util.Properties
import java.io.{File, FileInputStream}
import org.broadinstitute.sting.utils.text.XReadLines
import org.broadinstitute.sting.queue.{QException, QArguments}

/**
 * Syntactic sugar for filling in a pipeline using a Scala script.
 */
object Pipeline extends Logging
{
  var inputPaths = List.empty[File]
  // TODO: Stop using globals and wrap in an execution environment.  Will need when overriding values per command.
  var useBsub = false
  var dryRun = false
  private var argMap = Map.empty[String, String]
  private var rules = List.empty[QRule]

  /**
   * Sugar that allows addRule( inputs -> outputs, command )
   */
  def addRule(rule: (Any, Any), command: String): Unit = {
    addRule(rule._1, rule._2, command)
  }

  private def addRule(inputs: Any, outputs: Any, command: String): Unit = {
    rules :::= List(new QRule(getFiles(inputs), getFiles(outputs), new QCommand(command)))
  }

  def run(args: Array[String], inputs: Any, outputs: Any) = {
    QArguments.parseArgs(args)

    var inputFiles = getFiles(inputs)
    var outputFiles = getFiles(outputs)

    var grapher = new JobGrapher(inputPaths.map(_.getCanonicalPath), argMap, rules, inputFiles, outputFiles)

    val scheduler = useBsub match {
      case false => new SimpleJobScheduler(grapher.jobs)
      case true => new DispatchJobScheduler(grapher.jobs)
    }

    scheduler.runJobs
  }

  /**
   * Parses files passed in various sugar forms into a List[QFile]
   */
  private def getFiles(files: Any) : List[QFile] = {
    files match {
      case null => List.empty[QFile]
      case Nil => List.empty[QFile]
      case path: String => List(new QFile(path))
      case file: QFile => List(file)
      // Any List or Tuple add the members to this list
      case product: Product => {
        var list = List.empty[QFile]
        for (fileList <- product.productIterator.toList.map(getFiles(_))) {
          list :::= fileList
        }
        list
      }
      case x => throw new QException("Unknown file type: " + x)
    }
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
      for (line <- new XReadLines(file).iterator) {
        addFile(line)
      }
    } else {
      inputPaths = inputPaths ::: List(file)
    }
  }
}
