package org.broadinstitute.sting.queue

import tools.nsc.MainGenericRunner
import org.broadinstitute.sting.queue.util.ClasspathUtils
import collection.mutable.ListBuffer
import org.broadinstitute.sting.queue.util.Logging

object QCommandLine extends Application with Logging {
  var usage = """usage: java -jar Queue.jar [-P name=value] [-P file.properties] [-I input.file] [-I input_files.list] [-bsub] [-bsubWait] [-dry] [-debug] -S pipeline.scala"""

  override def main(args: Array[String]) = {
    val qArgs: QArguments = try {
      new QArguments(args)
    } catch {
      case exception => {
        println(exception)
        println(usage)
        System.exit(-1)
      }
      null
    }

    logger.debug("starting")

    if (qArgs.scripts.size == 0) {
      println("Error: Missing script")
      println(usage)
      System.exit(-1)
    }

    // NOTE: Something in MainGenericRunner is exiting the VM.
    if (qArgs.scripts.size != 1) {
      println("Error: Only one script can be run at a time")
      println(usage)
      System.exit(-1)
    }

    val newArgs = new ListBuffer[String]
    newArgs.appendAll(args)
    QArguments.strip(newArgs, "-S")
    newArgs.prepend("-nocompdaemon", "-classpath", ClasspathUtils.manifestAwareClassPath, qArgs.scripts.head)
    MainGenericRunner.main(newArgs.toArray)

    // NOTE: This line is not reached because the MainGenericRunner exits the VM.
    logger.debug("exiting")
  }
}
