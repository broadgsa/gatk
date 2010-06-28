package org.broadinstitute.sting.queue

import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.queue.engine.QGraph

/**
 * Syntactic sugar for filling in a pipeline using a Scala script.
 */
object QScript {
  // Type aliases so users don't have to import
  type File = java.io.File
  type Input = org.broadinstitute.sting.commandline.Input
  type Output = org.broadinstitute.sting.commandline.Output
  type CommandLineFunction = org.broadinstitute.sting.queue.function.CommandLineFunction
  type GatkFunction = org.broadinstitute.sting.queue.function.gatk.GatkFunction
  type ScatterGatherableFunction = org.broadinstitute.sting.queue.function.scattergather.ScatterGatherableFunction
  type Scatter = org.broadinstitute.sting.queue.function.scattergather.Scatter
  type Gather = org.broadinstitute.sting.queue.function.scattergather.Gather
  type BamGatherFunction = org.broadinstitute.sting.queue.function.scattergather.BamGatherFunction
  type SimpleTextGatherFunction = org.broadinstitute.sting.queue.function.scattergather.SimpleTextGatherFunction

  // The arguments for executing pipelines
  private var qArgs: QArguments = _

  // A default pipeline.  Can also use multiple 'new Pipeline()'
  private val pipeline = new Pipeline

  /**
   * Initializes the QArguments and returns a list of the rest of the user args.
   */
  def setArgs(params: Array[String]) = {
    qArgs = new QArguments(params)
    qArgs.userArgs
  }

  /**
   * Returns a list of files that were specified with "-I <file>" on the command line
   * or inside a .list file.
   */
  def inputs(extension: String) = qArgs.inputPaths.filter(_.getName.endsWith(extension))

  /**
   * Exchanges the extension on a file.
   */
  def swapExt(file: File, oldExtension: String, newExtension: String) =
    new File(file.getName.stripSuffix(oldExtension) + newExtension)

  /**
   * Adds one or more command line functions for dispatch later during run()
   */
  def add(functions: CommandLineFunction*) = pipeline.add(functions:_*)

  /**
   * Sets the @Input and @Output values for all the functions
   */
  def setParams(): Unit = pipeline.setParams()

  /**
   * Sets the @Input and @Output values for a single function
   */
  def setParams(function: CommandLineFunction): Unit = pipeline.setParams(function)

  /**
   * Executes functions that have been added to the pipeline.
   */
  def run() = pipeline.run()


  /**
   * Encapsulates a set of functions to run together.
   */
  protected class Pipeline {
    private var functions = List.empty[CommandLineFunction]

    /**
     * Adds one or more command line functions for dispatch later during run()
     */
    def add(functions: CommandLineFunction*) =
      this.functions :::= List(functions:_*)

    /**
     * Sets the @Input and @Output values for all the functions
     */
    def setParams(): Unit =
      for (function <- functions) setParams(function)

    /**
     * Sets the @Input and @Output values for a single function
     */
    def setParams(function: CommandLineFunction): Unit =
      function.properties = qArgs.properties

    /**
     * Executes functions that have been added to the pipeline.
     */
    def run() = {
      val qGraph = new QGraph
      qGraph.dryRun = qArgs.dryRun
      qGraph.bsubAllJobs = qArgs.bsubAllJobs
      qGraph.bsubWaitJobs = qArgs.bsubWaitJobs
      qGraph.properties = qArgs.properties
      for (function <- functions)
        qGraph.add(function)
      qGraph.fillIn
      qGraph.run
    }
  }
}
