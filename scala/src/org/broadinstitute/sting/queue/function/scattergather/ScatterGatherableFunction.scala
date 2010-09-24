package org.broadinstitute.sting.queue.function.scattergather

import java.io.File
import org.broadinstitute.sting.queue.util._
import org.broadinstitute.sting.commandline.ArgumentSource
import org.broadinstitute.sting.queue.function.CommandLineFunction
import com.rits.cloning.Cloner

/**
 * A function that can be run faster by splitting it up into pieces and then joining together the results.
 */
trait ScatterGatherableFunction extends CommandLineFunction {

  /** Number of parts to scatter the function into" */
  var scatterCount: Int = 1

  /** scatter gather directory */
  var scatterGatherDirectory: File = _

  /** cleanup temporary directories */
  var cleanupTempDirectories = false

  /** Class to use for creating temporary directories.  Defaults to CreateTempDirsFunction. */
  var createTempDirsClass: Class[_ <: CreateTempDirsFunction] = _

  /** Class to use for scattering.  Defaults to the annotation used in the @Scatter tag. */
  var scatterClass: Class[_ <: ScatterFunction] = _

  /**
   * Function that returns the class to use for gathering a directory.  If it returns null then @Gather annotation will be used.
   * @param gatherField Field that is to be gathered.
   * @return The class of the GatherFunction to be used or null.
   */
  var gatherClass: PartialFunction[ArgumentSource, Class[_ <: GatherFunction]] = _

  /** Class to use for removing temporary directories.  Defaults to CleanupTempDirsFunction. */
  var cleanupTempDirsClass: Class[_ <: CleanupTempDirsFunction] = _

  /**
   * Allows external modification of the CreateTempDirsFunction that will create the temporary directories.
   * @param initializeFunction The function that will create the temporary directories.
   * @param inputFields The input fields that the original function was dependent on.
   */
  var setupInitializeFunction: PartialFunction[(CreateTempDirsFunction, List[ArgumentSource]), Unit] = _

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
   * @param cloneFunction The clone of this ScatterGatherableFunction
   * @param index The one based index (from 1..scatterCount inclusive) of the scatter piece.
   */
  var setupCloneFunction: PartialFunction[(ScatterGatherableFunction, Int), Unit] = _

  /**
   * Allows external modification of the CleanupTempDirsFunction that will remove the temporary directories.
   * @param cleanupFunction The function that will remove the temporary directories.
   * @param gatherFunctions The functions that will gather up the original output fields.
   * @param outputFields The output fields that the original function was dependent on.
   */
  var setupCleanupFunction: PartialFunction[(CleanupTempDirsFunction, Map[ArgumentSource, GatherFunction], List[ArgumentSource]), Unit] = _

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
   * @return List[CommandLineFunction] to run instead of this function.
   */
  def generateFunctions() = {
    var functions = List.empty[CommandLineFunction]
    var tempDirectories = List.empty[File]

    // Only depend on input fields that have a value
    val inputFieldsWithValues = this.inputFields.filter(hasFieldValue(_))
    // Only gather up fields that will have a value
    val outputFieldsWithValues = this.outputFields.filter(hasFieldValue(_))

    // Create the scatter function based on @Scatter
    val scatterFunction = this.newScatterFunction(this.scatterField)
    initScatterFunction(scatterFunction, this.scatterField)
    tempDirectories :+= scatterFunction.commandDirectory
    functions :+= scatterFunction

    // Create the gather functions for each output field
    var gatherFunctions = Map.empty[ArgumentSource, GatherFunction]
    for (gatherField <- outputFieldsWithValues) {
      val gatherFunction = this.newGatherFunction(gatherField)
      initGatherFunction(gatherFunction, gatherField)
      tempDirectories :+= gatherFunction.commandDirectory
      functions :+= gatherFunction
      gatherFunctions += gatherField -> gatherFunction
    }

    // Create the clone functions for running the parallel jobs
    var cloneFunctions = List.empty[CommandLineFunction]
    for (i <- 1 to this.scatterCount) {
      val cloneFunction = this.newCloneFunction()
      initCloneFunction(cloneFunction, i)
      cloneFunctions :+= cloneFunction
      tempDirectories :+= cloneFunction.commandDirectory

      bindCloneFunctionScatter(scatterFunction, this.scatterField, cloneFunction, i)
      // For each each output field, change value to the scatterGatherTempDir dir and feed it into the gatherer
      for (gatherField <- outputFieldsWithValues)
        bindCloneFunctionGather(gatherFunctions(gatherField), gatherField, cloneFunction, i)
    }
    functions ++= cloneFunctions

    // Create a function to create all of the scatterGatherTempDir directories.
    // All of its inputs are the inputs of the original function.
    val initializeFunction = this.newInitializeFunction()
    initInitializeFunction(initializeFunction, inputFieldsWithValues)

    // Create a function that will remove any temporary items
    // All of its inputs are the outputs of the original function.
    var cleanupFunction = newCleanupFunction()
    initCleanupFunction(cleanupFunction, gatherFunctions, outputFieldsWithValues)

    // Set the temporary directories, for the initialize function as outputs for scatter and cleanup as inputs.
    initializeFunction.tempDirectories = tempDirectories
    scatterFunction.tempDirectories = tempDirectories
    cleanupFunction.tempDirectories = tempDirectories

    functions +:= initializeFunction
    if (this.cleanupTempDirectories)
      functions :+= cleanupFunction

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
   * Creates a new initialize CreateTempDirsFunction that will create the temporary directories.
   * @return A CreateTempDirsFunction that will create the temporary directories.
   */
  protected def newInitializeFunction(): CreateTempDirsFunction = {
    if (createTempDirsClass != null)
      this.createTempDirsClass.newInstance
    else
      new CreateTempDirsFunction
  }

  /**
   * Initializes the CreateTempDirsFunction that will create the temporary directories.
   * The initializeFunction qSettings is set so that the CreateTempDirsFunction runs with the same prefix, etc. as this ScatterGatherableFunction.
   * The initializeFunction commandDirectory is set so that the function runs in the directory as this ScatterGatherableFunction.
   * The initializeFunction is modified to become dependent on the input files for this ScatterGatherableFunction.
   * Calls setupInitializeFunction with initializeFunction.
   * @param initializeFunction The function that will create the temporary directories.
   * @param inputFields The input fields that the original function was dependent on.
   */
  protected def initInitializeFunction(initializeFunction: CreateTempDirsFunction, inputFields: List[ArgumentSource]) = {
    initializeFunction.qSettings = this.qSettings
    initializeFunction.commandDirectory = this.commandDirectory
    for (inputField <- inputFields)
      initializeFunction.originalInputs ++= this.getFieldFiles(inputField)
    if (this.setupInitializeFunction != null)
      if (this.setupInitializeFunction.isDefinedAt(initializeFunction, inputFields))
        this.setupInitializeFunction(initializeFunction, inputFields)
  }

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
   * The scatterFunction qSettings is set so that the ScatterFunction runs with the same prefix, etc. as this ScatterGatherableFunction.
   * The scatterFunction commandDirectory is set so that the function runs from a temporary directory under the scatterDirectory.
   * The scatterFunction has it's originalInput set with the file to be scattered into scatterCount pieces.
   * Calls scatterFunction.setOriginalFunction with this ScatterGatherableFunction.
   * Calls setupScatterFunction with scatterFunction.
   * @param scatterFunction The function that will create the scatter pieces in the temporary directories.
   * @param scatterField The input field being scattered.
   */
  protected def initScatterFunction(scatterFunction: ScatterFunction, scatterField: ArgumentSource) = {
    scatterFunction.qSettings = this.qSettings
    scatterFunction.commandDirectory = this.scatterGatherTempDir("scatter-" + scatterField.field.getName)
    scatterFunction.originalInput = this.getFieldFile(scatterField)
    scatterFunction.setOriginalFunction(this, scatterField)
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
   * The gatherFunction qSettings is set so that the GatherFunction runs with the same prefix, etc. as this ScatterGatherableFunction.
   * The gatherFunction commandDirectory is set so that the function runs from a temporary directory under the scatterDirectory.
   * The gatherFunction has it's originalOutput set with the file to be gathered from the scatterCount pieces.
   * Calls the gatherFunction.setOriginalFunction with this ScatterGatherableFunction.
   * Calls setupGatherFunction with gatherFunction.
   * @param gatherFunction The function that will merge the gather pieces from the temporary directories.
   * @param gatherField The output field being gathered.
   */
  protected def initGatherFunction(gatherFunction: GatherFunction, gatherField: ArgumentSource) = {
    gatherFunction.qSettings = this.qSettings
    gatherFunction.commandDirectory = this.scatterGatherTempDir("gather-" + gatherField.field.getName)
    gatherFunction.originalOutput = this.getFieldFile(gatherField)
    gatherFunction.setOriginalFunction(this, gatherField)
    gatherFunction.analysisName = this.analysisName
    if (this.setupGatherFunction != null)
      if (this.setupGatherFunction.isDefinedAt(gatherFunction, gatherField))
        this.setupGatherFunction(gatherFunction, gatherField)
  }

  /**
   * Creates a new clone of this ScatterGatherableFunction, setting the scatterCount to 1 so it doesn't infinitely scatter.
   * @return A clone of this ScatterGatherableFunction
   */
  protected def newCloneFunction(): ScatterGatherableFunction = {
    val cloneFunction = ScatterGatherableFunction.cloner.deepClone(this)
    // Make sure clone doesn't get scattered
    cloneFunction.scatterCount = 1
    cloneFunction
  }

  /**
   * Initializes the cloned function created by newCloneFunction() by setting it's commandDirectory to a temporary directory under scatterDirectory.
   * Calls setupCloneFunction with cloneFunction.
   * @param cloneFunction The clone of this ScatterGatherableFunction
   * @param index The one based index (from 1..scatterCount inclusive) of the scatter piece.
   */
  protected def initCloneFunction(cloneFunction: ScatterGatherableFunction, index: Int) = {
    cloneFunction.commandDirectory = this.scatterGatherTempDir("temp-"+index)
    if (this.setupCloneFunction != null)
      if (this.setupCloneFunction.isDefinedAt(cloneFunction, index))
        this.setupCloneFunction(cloneFunction, index)
        cloneFunction.isGather = false
  }

  /**
   * Joins a piece of the ScatterFunction output to the cloned function's input.
   * The input of the clone is changed to be in the output directory of the clone.
   * The scatter function piece is added as an output of the scatterFunction.
   * The clone function's original input is changed to use the piece from the output directory.
   * Finally the scatterFunction.setCloneFunction is called with the clone of this ScatterGatherableFunction.
   * @param scatterFunction Function that will create the pieces including the piece that will go to cloneFunction.
   * @param scatterField The field to be scattered.
   * @param cloneFunction Clone of this ScatterGatherableFunction.
   * @param index The one based index (from 1..scatterCount inclusive) of the scatter piece.
   */
  protected def bindCloneFunctionScatter(scatterFunction: ScatterFunction, scatterField: ArgumentSource, cloneFunction: ScatterGatherableFunction, index: Int) = {
    // Reset the input of the clone to the the scatterGatherTempDir dir and add it as an output of the scatter
    val scatterPart = IOUtils.resetParent(cloneFunction.commandDirectory, scatterFunction.originalInput)
    scatterFunction.scatterParts :+= scatterPart
    cloneFunction.setFieldValue(scatterField, scatterPart)
    scatterFunction.setCloneFunction(cloneFunction, index, scatterField)
  }

  /**
   * Joins the cloned function's output as a piece of the GatherFunction's input.
   * Finally the scatterFunction.setCloneFunction is called with the clone of this ScatterGatherableFunction.
   * @param cloneFunction Clone of this ScatterGatherableFunction.
   * @param gatherFunction Function that will create the pieces including the piece that will go to cloneFunction.
   * @param gatherField The field to be gathered.
   */
  protected def bindCloneFunctionGather(gatherFunction: GatherFunction, gatherField: ArgumentSource, cloneFunction: ScatterGatherableFunction, index: Int) = {
    val gatherPart = cloneFunction.resetFieldFile(gatherField, cloneFunction.commandDirectory)
    gatherFunction.gatherParts :+= gatherPart
    gatherFunction.setCloneFunction(cloneFunction, index, gatherField)
  }

  /**
   * Creates a new function that will remove the temporary directories.
   * @return A CleanupTempDirs function that will remove the temporary directories.
   */
  protected def newCleanupFunction(): CleanupTempDirsFunction = {
    if (cleanupTempDirsClass != null)
      this.cleanupTempDirsClass.newInstance
    else
      new CleanupTempDirsFunction
  }

  /**
   * Initializes the CleanupTempDirsFunction created by newCleanupFunction() that will remove the temporary directories.
   * The cleanupFunction qSettings is set so that the CleanupTempDirsFunction runs with the same prefix, etc. as this ScatterGatherableFunction.
   * The cleanupFunction commandDirectory is set so that the function runs in the directory as this ScatterGatherableFunction.
   * The initializeFunction is modified to become dependent on the output files for this ScatterGatherableFunction.
   * Calls setupCleanupFunction with cleanupFunction.
   * @param cleanupFunction The function that will remove the temporary directories.
   * @param gatherFunctions The functions that will gather up the original output fields.
   * @param outputFields The output fields that the original function was dependent on.
   */
  protected def initCleanupFunction(cleanupFunction: CleanupTempDirsFunction, gatherFunctions: Map[ArgumentSource, GatherFunction], outputFields: List[ArgumentSource]) = {
    cleanupFunction.qSettings = this.qSettings
    cleanupFunction.commandDirectory = this.commandDirectory
    for (gatherField <- outputFields)
      cleanupFunction.originalOutputs += gatherFunctions(gatherField).originalOutput
    if (this.setupCleanupFunction != null)
      if (this.setupCleanupFunction.isDefinedAt(cleanupFunction, gatherFunctions, outputFields))
        this.setupCleanupFunction(cleanupFunction, gatherFunctions, outputFields)
  }

  /**
   * Returns a temporary directory under this scatter gather directory.
   * @param Sub directory under the scatter gather directory.
   * @return temporary directory under this scatter gather directory.
   */
  private def scatterGatherTempDir(subDir: String) = IOUtils.subDir(this.scatterGatherDirectory, this.jobName + "-" + subDir)
}

/**
 * A function that can be run faster by splitting it up into pieces and then joining together the results.
 */
object ScatterGatherableFunction {
  /** Used to deep clone a ScatterGatherableFunction. */
  private lazy val cloner = new Cloner
}
