package org.broadinstitute.sting.queue

import engine.Pipeline
import tools.nsc.MainGenericRunner
import org.broadinstitute.sting.queue.util.ClasspathUtils
import org.apache.log4j._
import collection.mutable.ListBuffer
import org.broadinstitute.sting.queue.util.Logging

object QCommandLine extends Application with Logging {
  var usage = """usage: java -jar Queue.jar [ -P name=value ] [ -P file.properties ] [ -I input.file ] [ -I input_files.list ] [ -bsub ] [ -dry ] [ -debug ] -S pipeline.scala"""

  override def main(args: Array[String]) = {
    try {
      QArguments.parseArgs(args.clone)
    } catch {
      case e: Exception => {
        println(usage)
        System.exit(-1)
      }
    }

    var newArgs = new ListBuffer[String]
    newArgs.appendAll(args)

    QArguments.strip(newArgs, "-S")

    logger.debug("starting")

    if (QArguments.scripts.size == 0) {
      println("Error: Missing script")
      println(usage)
      System.exit(-1)
    }

    if (Pipeline.inputPaths.size == 0) {
      println("Error: No inputs specified")
      println(usage)
      System.exit(-1)
    }

    for (script <- QArguments.scripts) {
      var clone = newArgs.clone
      clone.prepend("-nocompdaemon", "-classpath", ClasspathUtils.manifestAwareClassPath, script)
      MainGenericRunner.main(clone.toArray)
    }

    logger.debug("exiting")
  }
}
