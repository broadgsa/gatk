/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.queue.function.scattergather

import java.io.File
import org.broadinstitute.gatk.queue.util._
import org.broadinstitute.gatk.utils.commandline.{Gatherer, Gather, ArgumentSource}
import org.broadinstitute.gatk.queue.function.{QFunction, CommandLineFunction}
import org.broadinstitute.gatk.queue.QException
import org.broadinstitute.gatk.utils.io.IOUtils
import collection.immutable.ListMap

/**
 * A function that can be run faster by splitting it up into pieces and then joining together the results.
 */
trait ScatterGatherableFunction extends CommandLineFunction {

  /** Maximum number of parts to scatter the function into. */
  var scatterCount: Int = 1

  /** scatter gather directory */
  var scatterGatherDirectory: File = _

  /** Class to use for scattering.  Defaults to the annotation used in the @Scatter tag. */
  var scatterClass: Class[_ <: ScatterFunction] = _

  /**
   * Function that returns the class to use for gathering a directory.  If it returns null then @Gather annotation will be used.
   * PartialFunction param gatherField Field that is to be gathered.
   * @return The class of the GatherFunction to be used or null.
   */
  var gatherClass: PartialFunction[ArgumentSource, Class[_ <: GatherFunction]] = _

  /**
   * Allows external modification of the ScatterFunction that will create the scatter pieces in the temporary directories.
   * PartialFunction param scatterFunction The function that will create the scatter pieces in the temporary directories.
   */
  var setupScatterFunction: PartialFunction[ScatterFunction, Unit] = _

  /**
   * Allows external modification of the GatherFunction that will collect the gather pieces in the temporary directories.
   * PartialFunction param gatherFunction The function that will merge the gather pieces from the temporary directories.
   * PartialFunction param gatherField The output field being gathered.
   */
  var setupGatherFunction: PartialFunction[(GatherFunction, ArgumentSource), Unit] = _

  /**
   * Allows external modification of the cloned function.
   * PartialFunction param cloneFunction A clone wrapper of this ScatterGatherableFunction
   * PartialFunction param index The one based index (from 1..scatterCount inclusive) of the scatter piece.
   */
  var setupCloneFunction: PartialFunction[(CloneFunction, Int), Unit] = _

  /**
   * Returns true if the function is ready to be scatter / gathered.
   * The base implementation checks if the scatter count is greater than one,
   * and that the scatter function can scatter this instance.
   * @return true if the function is ready to be scatter / gathered.
   */
  def scatterGatherable = this.scatterCount > 1 && scatterFunction.isScatterGatherable

  /**
   * Sets the scatter gather directory to the command directory if it is not already set.
   */
  override def freezeFieldValues() {
    super.freezeFieldValues()

    if (this.scatterGatherDirectory == null) {
      if (qSettings.jobScatterGatherDirectory != null) {
        this.scatterGatherDirectory = IOUtils.absolute(qSettings.jobScatterGatherDirectory)
      } else {
        this.scatterGatherDirectory = IOUtils.absolute(this.commandDirectory, ".queue/scatterGather")
      }
    }
  }

  /**
   * The scatter function.
   */
  private lazy val scatterFunction = {
    // Only depend on input fields that have a value
    val inputFieldsWithValues = this.inputFields.filter(hasFieldValue(_))
    val inputFiles = inputFieldsWithValues.flatMap(getFieldFiles(_)).toSet

    val scatterFunction = newScatterFunction()
    this.copySettingsTo(scatterFunction)
    scatterFunction.originalFunction = this
    scatterFunction.originalInputs = inputFiles
    scatterFunction.commandDirectory = this.scatterGatherCommandDir("scatter")
    scatterFunction.jobOutputFile = new File("scatter.out")
    scatterFunction.addOrder = this.addOrder :+ 1
    scatterFunction.isIntermediate = true

    initScatterFunction(scatterFunction)
    scatterFunction.absoluteCommandDirectory()
    scatterFunction.init()
    scatterFunction
  }

  /**
   * Returns a list of scatter / gather and clones of this function
   * that can be run in parallel to produce the same output as this
   * command line function.
   * @return Seq[QFunction] to run instead of this function.
   */
  def generateFunctions() = {
    // Ask the scatter function how many clones to create.
    val numClones = scatterFunction.scatterCount

    // Create the gather functions for each output field
    var gatherFunctions = ListMap.empty[ArgumentSource, GatherFunction]
    var gatherOutputs = ListMap.empty[ArgumentSource, File]
    var gatherAddOrder = numClones + 2

    // Only track fields that will have an output file
    val outputFieldsWithValues = this.outputFields.
      filter(hasFieldValue(_)).
      filter(gatherField => !IOUtils.isSpecialFile(getFieldFile(gatherField)))

    for (gatherField <- outputFieldsWithValues) {
      gatherOutputs += gatherField -> getFieldFile(gatherField)
    }

    // Only gather fields that are @Gather(enabled=true)
    val outputFieldsWithGathers = outputFieldsWithValues.filter(hasGatherFunction(_))

    for (gatherField <- outputFieldsWithGathers) {
      val gatherOutput = gatherOutputs(gatherField)

      val gatherFunction = this.newGatherFunction(gatherField)
      this.copySettingsTo(gatherFunction)
      gatherFunction.originalFunction = this
      gatherFunction.originalOutput = gatherOutput
      gatherFunction.commandDirectory = this.scatterGatherCommandDir("gather-" + gatherField.field.getName)
      gatherFunction.jobOutputFile = new File("gather-" + gatherOutput.getName + ".out")
      gatherFunction.addOrder = this.addOrder :+ gatherAddOrder

      initGatherFunction(gatherFunction, gatherField)
      gatherFunction.absoluteCommandDirectory()
      gatherFunction.init()

      gatherFunctions += gatherField -> gatherFunction

      gatherAddOrder += 1
    }

    // Create the clone functions for running the parallel jobs
    var cloneFunctions = Seq.empty[CloneFunction]
    val dirFormat = "temp_%%0%dd_of_%d".format(numClones.toString.length(), numClones)
    for (i <- 1 to numClones) {
      val cloneFunction = this.newCloneFunction()

      this.copySettingsTo(cloneFunction)
      cloneFunction.originalFunction = this
      cloneFunction.analysisName = this.analysisName
      cloneFunction.cloneIndex = i
      cloneFunction.cloneCount = numClones
      cloneFunction.commandDirectory = this.scatterGatherCommandDir(dirFormat.format(i))
      cloneFunction.jobOutputFile = if (IOUtils.isSpecialFile(this.jobOutputFile)) this.jobOutputFile else new File(this.jobOutputFile.getName)
      if (this.jobErrorFile != null)
        cloneFunction.jobErrorFile = if (IOUtils.isSpecialFile(this.jobErrorFile)) this.jobErrorFile else new File(this.jobErrorFile.getName)
      // jic the "local" dir is actually on the network, create different sub local directories for each clone.
      // This might be better handled with a hook that allows clones to create unique file names. Right now no hook
      // like freezeFieldValues exists for specifying per cloneFunction fields.
      cloneFunction.jobLocalDir = this.scatterGatherLocalDir(dirFormat.format(i))
      cloneFunction.addOrder = this.addOrder :+ (i+1)
      cloneFunction.isIntermediate = true

      // Setup the fields on the clone function, outputting each as a relative file in the sg directory.
      scatterFunction.initCloneInputs(cloneFunction, i)
      for (gatherField <- outputFieldsWithValues) {
        val gatherPart = new File(gatherOutputs(gatherField).getName)
        cloneFunction.setFieldValue(gatherField, gatherPart)
      }

      // Allow the script writer to change the paths to the files.
      initCloneFunction(cloneFunction, i)

      // If the command directory is relative, insert the run directory ahead of it.
      cloneFunction.absoluteCommandDirectory()

      // Allow the scatter function to set the specific input for this clone
      scatterFunction.bindCloneInputs(cloneFunction, i)

      // Set each of the clone outputs to be absolute paths.
      for (gatherField <- outputFieldsWithValues) {
        val gatherPart = IOUtils.absolute(cloneFunction.commandDirectory, cloneFunction.getFieldFile(gatherField))
        cloneFunction.setFieldValue(gatherField, gatherPart)
      }

      // For the outputs that are being gathered add this clone's output to be gathered.
      for (gatherField <- outputFieldsWithGathers) {
        gatherFunctions(gatherField).gatherParts :+= cloneFunction.getFieldFile(gatherField)
      }

      cloneFunctions :+= cloneFunction
    }

    // Track the functions starting with the scatter function.
    var functions: Seq[QFunction] = Seq(scatterFunction) ++ cloneFunctions ++ gatherFunctions.values

    // Make all log file paths absolute.
    for (function <- functions) {
      function.jobOutputFile = IOUtils.absolute(function.commandDirectory, function.jobOutputFile)
      if (function.jobErrorFile != null)
        function.jobErrorFile = IOUtils.absolute(function.commandDirectory, function.jobErrorFile)
    }

    val jobOutputGather = gatherLogFile(_.jobOutputFile, functions, gatherAddOrder)
    if (this.jobErrorFile != null) {
      val jobErrorGather = gatherLogFile(_.jobErrorFile, functions, gatherAddOrder + 1)
      functions :+= jobErrorGather
    }
    functions :+= jobOutputGather

    // Return all the various created functions.
    functions
  }

  /**
   * Creates a new ScatterFunction.
   * @return A ScatterFunction instantiated scatterClass
   */
  protected def newScatterFunction(): ScatterFunction = {
    if (this.scatterClass == null)
      throw new QException("scatterClass is null.")
    this.scatterClass.newInstance.asInstanceOf[ScatterFunction]
  }

  /**
   * Initializes the ScatterFunction created by newScatterFunction() that will create the scatter pieces in the temporary directories.
   * Calls setupScatterFunction with scatterFunction.
   * @param scatterFunction The function that will create the scatter pieces in the temporary directories.
   */
  protected def initScatterFunction(scatterFunction: ScatterFunction) {
    if (this.setupScatterFunction != null)
      if (this.setupScatterFunction.isDefinedAt(scatterFunction))
        this.setupScatterFunction(scatterFunction)
  }

  /**
   * Returns true if the field should be gathered.
   * @param gatherField Field that defined @Gather.
   * @return true if the field should be gathered.
   */
  protected def hasGatherFunction(gatherField: ArgumentSource) : Boolean = {
    // Check if there is a function that will return the gather class for this field.
    if (this.gatherClass != null && this.gatherClass.isDefinedAt(gatherField))
        true

    // Check for an annotation defining the gather class.
    else if (ReflectionUtils.hasAnnotation(gatherField.field, classOf[Gather]))
      ReflectionUtils.getAnnotation(gatherField.field, classOf[Gather]).enabled

    // Nothing else to disable this field.
    else
      true
  }

  /**
   * Creates a new GatherFunction for the gatherField.
   * @param gatherField Field that defined @Gather.
   * @return A GatherFunction instantiated from @Gather.
   */
  protected def newGatherFunction(gatherField: ArgumentSource) : GatherFunction = {
    var gatherClass: Class[_] = null

    // Check if there is a function that will return the gather class for this field.
    if (this.gatherClass != null)
      if (this.gatherClass.isDefinedAt(gatherField))
        gatherClass = this.gatherClass(gatherField)

    // Check for an annotation defining the gather class.
    if (gatherClass == null) {
      if (ReflectionUtils.hasAnnotation(gatherField.field, classOf[Gather])) {
        gatherClass = ReflectionUtils.getAnnotation(gatherField.field, classOf[Gather]).value
      } else {
        throw new QException("Missing @Gather annotation on %s".format(gatherField.field))
      }
    }

    if (gatherClass == classOf[GatherFunction]) {
      throw new QException("@Gather did not specify class type on %s".format(gatherField.field))
    } else if (classOf[GatherFunction].isAssignableFrom(gatherClass)) {
      gatherClass.newInstance.asInstanceOf[GatherFunction]
    } else if (classOf[Gatherer].isAssignableFrom(gatherClass)) {
      new GathererFunction(gatherClass.asSubclass(classOf[Gatherer]))
    } else {
      throw new QException("Unsupported @Gather class type on %s: %s".format(gatherField.field, gatherClass))
    }
  }

  /**
   * Initializes the GatherFunction created by newGatherFunction() that will collect the gather pieces in the temporary directories.
   * Calls the gatherFunction.setOriginalFunction with this ScatterGatherableFunction.
   * Calls setupGatherFunction with gatherFunction.
   * @param gatherFunction The function that will merge the gather pieces from the temporary directories.
   * @param gatherField The output field being gathered.
   */
  protected def initGatherFunction(gatherFunction: GatherFunction, gatherField: ArgumentSource) {
    if (this.setupGatherFunction != null)
      if (this.setupGatherFunction.isDefinedAt(gatherFunction, gatherField))
        this.setupGatherFunction(gatherFunction, gatherField)
  }

  /**
   * Creates a new clone of this ScatterGatherableFunction, setting the scatterCount to 1 so it doesn't infinitely scatter.
   * @return An uninitialized clone wrapper for ScatterGatherableFunction
   */
  protected def newCloneFunction() = new CloneFunction

  /**
   * Calls setupCloneFunction with cloneFunction.
   * @param cloneFunction The clone of this ScatterGatherableFunction
   * @param index The one based index (from 1..scatterCount inclusive) of the scatter piece.
   */
  protected def initCloneFunction(cloneFunction: CloneFunction, index: Int) {
    if (this.setupCloneFunction != null)
      if (this.setupCloneFunction.isDefinedAt(cloneFunction, index))
        this.setupCloneFunction(cloneFunction, index)
  }

  /**
   * Gathers up the logs files from other functions.
   * @param logFile Takes the QFunction and return the log file.
   * @param functions The functions for which the logs will be concatenated.
   * @param addOrder The order this function should be added in the graph.
   */
  private def gatherLogFile(logFile: (QFunction) => File, functions: Seq[QFunction], addOrder: Int) = {
    val gatherLogFunction = new ConcatenateLogsFunction
    this.copySettingsTo(gatherLogFunction)
    gatherLogFunction.logs = functions.map(logFile).filter(_ != null)
    gatherLogFunction.jobOutputFile = logFile(this)
    gatherLogFunction.commandDirectory = this.scatterGatherCommandDir()
    gatherLogFunction.addOrder = this.addOrder :+ addOrder
    gatherLogFunction.isIntermediate = false
    gatherLogFunction
  }

  /**
   * Returns a temporary directory under this scatter gather directory.
   * @param subDir directory under the scatter gather directory.
   * @return temporary directory under this scatter gather directory.
   */
  private def scatterGatherCommandDir(subDir: String = "") = IOUtils.absolute(this.scatterGatherDirectory, this.jobName + "-sg/" + subDir)

  /**
   * Returns a sub directory under this job local directory.
   * @param subDir directory under the job local directory.
   * @return absolute path to a directory under the original job local directory.
   */
  private def scatterGatherLocalDir(subDir: String = "") = IOUtils.absolute(this.jobLocalDir, this.jobName + "-sg/" + subDir)
}
