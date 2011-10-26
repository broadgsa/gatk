/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.queue.function.scattergather

import java.io.File
import org.broadinstitute.sting.queue.util._
import org.broadinstitute.sting.commandline.{Gatherer, Gather, ArgumentSource}
import org.broadinstitute.sting.queue.function.{QFunction, CommandLineFunction}
import org.broadinstitute.sting.queue.QException
import org.broadinstitute.sting.utils.io.IOUtils

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
   * @param gatherField Field that is to be gathered.
   * @return The class of the GatherFunction to be used or null.
   */
  var gatherClass: PartialFunction[ArgumentSource, Class[_ <: GatherFunction]] = _

  /**
   * Allows external modification of the ScatterFunction that will create the scatter pieces in the temporary directories.
   * @param scatterFunction The function that will create the scatter pieces in the temporary directories.
   */
  var setupScatterFunction: PartialFunction[ScatterFunction, Unit] = _

  /**
   * Allows external modification of the GatherFunction that will collect the gather pieces in the temporary directories.
   * @param gatherFunction The function that will merge the gather pieces from the temporary directories.
   * @param gatherField The output field being gathered.
   */
  var setupGatherFunction: PartialFunction[(GatherFunction, ArgumentSource), Unit] = _

  /**
   * Allows external modification of the cloned function.
   * @param cloneFunction A clone wrapper of this ScatterGatherableFunction
   * @param index The one based index (from 1..scatterCount inclusive) of the scatter piece.
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
        this.scatterGatherDirectory = IOUtils.absolute(this.commandDirectory, "queueScatterGather")
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
    scatterFunction.commandDirectory = this.scatterGatherTempDir("scatter")
    scatterFunction.isIntermediate = true
    scatterFunction.addOrder = this.addOrder :+ 1

    initScatterFunction(scatterFunction)
    scatterFunction.absoluteCommandDirectory()
    scatterFunction.init()
    scatterFunction
  }

  /**
   * Returns a list of scatter / gather and clones of this function
   * that can be run in parallel to produce the same output as this
   * command line function.
   * @return List[QFunction] to run instead of this function.
   */
  def generateFunctions() = {
    var functions = List.empty[QFunction]

    // Only gather up fields that will have a value
    val outputFieldsWithValues = this.outputFields.filter(hasFieldValue(_))

    // Create the scatter function based on @Scatter
    functions :+= scatterFunction

    // Ask the scatter function how many clones to create.
    val numClones = scatterFunction.scatterCount

    // List of the log files that are output by this function.
    var logFiles = List(this.jobOutputFile)
    if (this.jobErrorFile != null)
      logFiles :+= this.jobErrorFile

    // Create the gather functions for each output field
    var gatherFunctions = Map.empty[ArgumentSource, GatherFunction]
    var gatherOutputs = Map.empty[ArgumentSource, File]
    var gatherAddOrder = numClones + 2
    for (gatherField <- outputFieldsWithValues) {
      val gatherOutput = getFieldFile(gatherField)

      val gatherFunction = this.newGatherFunction(gatherField)
      this.copySettingsTo(gatherFunction)
      gatherFunction.originalFunction = this
      gatherFunction.originalOutput = gatherOutput
      gatherFunction.commandDirectory = this.scatterGatherTempDir("gather-" + gatherField.field.getName)
      // If this is a gather for a log file, make the gather intermediate just in case the log file name changes
      // Otherwise have the regular output function wait on the log files to gather
      if (isLogFile(gatherOutput)) {
        gatherFunction.isIntermediate = true
        // Only delete the log files if the original function is an intermediate
        // and the intermediate files are supposed to be deleted
        gatherFunction.deleteIntermediateOutputs = this.isIntermediate && this.deleteIntermediateOutputs
      } else {
        gatherFunction.originalLogFiles = logFiles
      }
      gatherFunction.addOrder = this.addOrder :+ gatherAddOrder

      initGatherFunction(gatherFunction, gatherField)
      gatherFunction.absoluteCommandDirectory()
      gatherFunction.init()

      functions :+= gatherFunction
      gatherFunctions += gatherField -> gatherFunction
      gatherOutputs += gatherField -> gatherOutput

      gatherAddOrder += 1
    }

    // Create the clone functions for running the parallel jobs
    var cloneFunctions = List.empty[CloneFunction]
    for (i <- 1 to numClones) {
      val cloneFunction = this.newCloneFunction()

      this.copySettingsTo(cloneFunction)
      cloneFunction.originalFunction = this
      cloneFunction.cloneIndex = i
      cloneFunction.commandDirectory = this.scatterGatherTempDir("temp-"+i)
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

      // Get absolute paths to the files and bind the sg functions to the clone function via the absolute paths.
      scatterFunction.bindCloneInputs(cloneFunction, i)
      for (gatherField <- outputFieldsWithValues) {
        val gatherPart = IOUtils.absolute(cloneFunction.commandDirectory, cloneFunction.getFieldFile(gatherField))
        cloneFunction.setFieldValue(gatherField, gatherPart)
        gatherFunctions(gatherField).gatherParts :+= gatherPart
      }

      cloneFunctions :+= cloneFunction
    }
    functions ++= cloneFunctions

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
        throw new QException("Missing @Gather annotation: " + gatherField.field)
      }
    }

    if (classOf[GatherFunction].isAssignableFrom(gatherClass)) {
      gatherClass.newInstance.asInstanceOf[GatherFunction]
    } else if (classOf[Gatherer].isAssignableFrom(gatherClass)) {
      new GathererFunction(gatherClass.asSubclass(classOf[Gatherer]))
    } else {
      throw new QException("Unsupported @Gather class type: " + gatherClass)
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
   * Returns a temporary directory under this scatter gather directory.
   * @param Sub directory under the scatter gather directory.
   * @return temporary directory under this scatter gather directory.
   */
  private def scatterGatherTempDir(subDir: String) = IOUtils.absolute(this.scatterGatherDirectory, this.jobName + "-sg/" + subDir)
}
