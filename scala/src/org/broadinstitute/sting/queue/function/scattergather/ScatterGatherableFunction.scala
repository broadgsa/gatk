package org.broadinstitute.sting.queue.function.scattergather

import java.io.File
import org.broadinstitute.sting.queue.util._
import org.broadinstitute.sting.commandline.ArgumentSource
import org.broadinstitute.sting.queue.function.{QFunction, CommandLineFunction}

/**
 * A function that can be run faster by splitting it up into pieces and then joining together the results.
 */
trait ScatterGatherableFunction extends CommandLineFunction {

  /** Number of parts to scatter the function into" */
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
   * @param scatterField The input field being scattered.
   */
  var setupScatterFunction: PartialFunction[(ScatterFunction, ArgumentSource), Unit] = _

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
   * and that the scatter field has a value.
   * @return true if the function is ready to be scatter / gathered.
   */
  def scatterGatherable = this.scatterCount > 1 && hasFieldValue(this.scatterField)

  /**
   * Returns a list of scatter / gather and clones of this function
   * that can be run in parallel to produce the same output as this
   * command line function.
   * @return List[QFunction] to run instead of this function.
   */
  def generateFunctions() = {
    var functions = List.empty[QFunction]

    // Only depend on input fields that have a value
    val inputFieldsWithValues = this.inputFields.filter(hasFieldValue(_))
    // Only gather up fields that will have a value
    val outputFieldsWithValues = this.outputFields.filter(hasFieldValue(_))
    // The field containing the file to split
    val originalInput = getFieldFile(scatterField)

    // Create the scatter function based on @Scatter
    val scatterFunction = this.newScatterFunction(this.scatterField)
    syncFunction(scatterFunction)
    scatterFunction.addOrder = this.addOrder :+ 1
    scatterFunction.commandDirectory = this.scatterGatherTempDir("scatter-" + scatterField.field.getName)
    scatterFunction.originalInput = originalInput
    scatterFunction.setOriginalFunction(this, scatterField)
    initScatterFunction(scatterFunction, this.scatterField)
    functions :+= scatterFunction

    // Create the gather functions for each output field
    var gatherFunctions = Map.empty[ArgumentSource, GatherFunction]
    var gatherOutputs = Map.empty[ArgumentSource, File]
    var gatherAddOrder = this.scatterCount + 2
    for (gatherField <- outputFieldsWithValues) {
      val gatherFunction = this.newGatherFunction(gatherField)
      val gatherOutput = getFieldFile(gatherField)
      syncFunction(gatherFunction)
      gatherFunction.addOrder = this.addOrder :+ gatherAddOrder
      gatherFunction.commandDirectory = this.scatterGatherTempDir("gather-" + gatherField.field.getName)
      gatherFunction.originalOutput = this.getFieldFile(gatherField)
      gatherFunction.setOriginalFunction(this, gatherField)
      initGatherFunction(gatherFunction, gatherField)
      functions :+= gatherFunction
      gatherFunctions += gatherField -> gatherFunction
      gatherOutputs += gatherField -> gatherOutput
      gatherAddOrder += 1
    }

    // Create the clone functions for running the parallel jobs
    var cloneFunctions = List.empty[CloneFunction]
    for (i <- 1 to this.scatterCount) {
      val cloneFunction = this.newCloneFunction()

      syncFunction(cloneFunction)
      cloneFunction.originalFunction = this
      cloneFunction.index = i
      cloneFunction.addOrder = this.addOrder :+ (i+1)
      cloneFunction.memoryLimit = this.memoryLimit

      // Setup the fields on the clone function, outputting each as a relative file in the sg directory.
      cloneFunction.commandDirectory = this.scatterGatherTempDir("temp-"+i)
      var scatterPart = new File(originalInput.getName)
      cloneFunction.setFieldValue(scatterField, scatterPart)
      for (gatherField <- outputFieldsWithValues) {
        val gatherPart = new File(gatherOutputs(gatherField).getName)
        cloneFunction.setFieldValue(gatherField, gatherPart)
      }

      // Allow the script writer to change the paths to the files.
      initCloneFunction(cloneFunction, i)

      // Get absolute paths to the files and bind the sg functions to the clone function via the absolute paths.
      scatterPart = IOUtils.subDir(cloneFunction.commandDirectory, cloneFunction.getFieldFile(scatterField))
      cloneFunction.setFieldValue(scatterField, scatterPart)
      scatterFunction.scatterParts :+= scatterPart
      for (gatherField <- outputFieldsWithValues) {
        val gatherPart = IOUtils.subDir(cloneFunction.commandDirectory, cloneFunction.getFieldFile(gatherField))
        cloneFunction.setFieldValue(gatherField, gatherPart)
        gatherFunctions(gatherField).gatherParts :+= gatherPart
      }

      cloneFunctions :+= cloneFunction
    }
    functions ++= cloneFunctions

    // Return all the various functions we created
    functions
  }

  /**
   * Sets the scatter gather directory to the command directory if it is not already set.
   */
  override def freezeFieldValues = {
    super.freezeFieldValues

    if (this.scatterGatherDirectory == null) {
      this.scatterGatherDirectory = qSettings.jobScatterGatherDirectory
      if (this.scatterGatherDirectory == null)
        this.scatterGatherDirectory = this.commandDirectory
    }
  }

  /**
   * Retrieves the scatter field from the first field that has the annotation @Scatter.
   */
  protected lazy val scatterField =
    this.inputFields.find(field => ReflectionUtils.hasAnnotation(field.field, classOf[Scatter])).get

  /**
   * Creates a new ScatterFunction for the scatterField.
   * @param scatterField Field that defined @Scatter.
   * @return A ScatterFunction instantiated from @Scatter or scatterClass if scatterClass was set on this ScatterGatherableFunction.
   */
  protected def newScatterFunction(scatterField: ArgumentSource): ScatterFunction = {
    var scatterClass = this.scatterClass
    if (scatterClass == null)
      scatterClass = ReflectionUtils.getAnnotation(scatterField.field, classOf[Scatter])
              .value.asSubclass(classOf[ScatterFunction])
    scatterClass.newInstance.asInstanceOf[ScatterFunction]
  }

  /**
   * Initializes the ScatterFunction created by newScatterFunction() that will create the scatter pieces in the temporary directories.
   * Calls setupScatterFunction with scatterFunction.
   * @param scatterFunction The function that will create the scatter pieces in the temporary directories.
   * @param scatterField The input field being scattered.
   */
  protected def initScatterFunction(scatterFunction: ScatterFunction, scatterField: ArgumentSource) = {
    if (this.setupScatterFunction != null)
      if (this.setupScatterFunction.isDefinedAt(scatterFunction, scatterField))
        this.setupScatterFunction(scatterFunction, scatterField)
  }

  /**
   * Creates a new GatherFunction for the gatherField.
   * @param gatherField Field that defined @Gather.
   * @return A GatherFunction instantiated from @Gather.
   */
  protected def newGatherFunction(gatherField: ArgumentSource) : GatherFunction = {
    var gatherClass: Class[_ <: GatherFunction] = null
    if (this.gatherClass != null)
      if (this.gatherClass.isDefinedAt(gatherField))
        gatherClass = this.gatherClass(gatherField)
    if (gatherClass == null)
      gatherClass = ReflectionUtils.getAnnotation(gatherField.field, classOf[Gather])
              .value.asSubclass(classOf[GatherFunction])
    gatherClass.newInstance.asInstanceOf[GatherFunction]
  }

  /**
   * Initializes the GatherFunction created by newGatherFunction() that will collect the gather pieces in the temporary directories.
   * Calls the gatherFunction.setOriginalFunction with this ScatterGatherableFunction.
   * Calls setupGatherFunction with gatherFunction.
   * @param gatherFunction The function that will merge the gather pieces from the temporary directories.
   * @param gatherField The output field being gathered.
   */
  protected def initGatherFunction(gatherFunction: GatherFunction, gatherField: ArgumentSource) = {
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
  protected def initCloneFunction(cloneFunction: CloneFunction, index: Int) = {
    if (this.setupCloneFunction != null)
      if (this.setupCloneFunction.isDefinedAt(cloneFunction, index))
        this.setupCloneFunction(cloneFunction, index)
  }

  /**
   * Copies standard values from this function to the just created function.
   * @param newFunction newly created function.
   */
  protected def syncFunction(newFunction: QFunction) = {
    newFunction.isIntermediate = this.isIntermediate
    newFunction.analysisName = this.analysisName
    newFunction.qSettings = this.qSettings
    newFunction.jobTempDir = this.jobTempDir
    newFunction.jobName = this.jobName
    newFunction match {
      case newCLFFunction: CommandLineFunction =>
        newCLFFunction.jobQueue = this.jobQueue
        newCLFFunction.jobProject = this.jobProject
      case _ => /* ignore */
    }
  }

  /**
   * Returns a temporary directory under this scatter gather directory.
   * @param Sub directory under the scatter gather directory.
   * @return temporary directory under this scatter gather directory.
   */
  private def scatterGatherTempDir(subDir: String) = IOUtils.subDir(this.scatterGatherDirectory, this.jobName + "-sg/" + subDir)
}
