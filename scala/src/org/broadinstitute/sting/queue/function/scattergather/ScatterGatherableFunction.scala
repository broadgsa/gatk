package org.broadinstitute.sting.queue.function.scattergather

import org.broadinstitute.sting.queue.function.CommandLineFunction
import java.lang.reflect.Field
import java.io.File
import org.broadinstitute.sting.queue.util._

trait ScatterGatherableFunction extends CommandLineFunction {
  @Internal
  var scatterCount: Int = 1

  def scatterField = this.inputFields.find(field => ReflectionUtils.hasAnnotation(field, classOf[Scatter])).get

  def scatterGatherable = {
    if (scatterCount < 2)
      false
    else if (!hasFieldValue(scatterField))
      false
    else
      true
  }

  def generateFunctions() = ScatterGatherableFunction.generateFunctions(this)
}

object ScatterGatherableFunction {
  private def generateFunctions(originalFunction: ScatterGatherableFunction) = {
    var functions = List.empty[CommandLineFunction]
    var tempDirectories = List.empty[File]

    // Create a function that will remove any temporary items
    var cleanupFunction = new CleanupTempDirsFunction
    cleanupFunction.jobNamePrefix = originalFunction.jobNamePrefix
    cleanupFunction.commandDirectory = originalFunction.commandDirectory

    // Find the field with @Scatter and its value
    var scatterField = originalFunction.scatterField
    val originalValue = originalFunction.getFieldValue(scatterField)

    // Create the scatter function based on @Scatter
    val scatterFunction = getScatterFunction(scatterField)
    scatterFunction.setOriginalFunction(originalFunction)
    scatterFunction.jobNamePrefix = originalFunction.jobNamePrefix
    scatterFunction.commandDirectory = originalFunction.temp("scatter-" + scatterField.getName)
    scatterFunction.originalInput = originalValue.asInstanceOf[scatterFunction.ScatterType]
    tempDirectories :+= scatterFunction.commandDirectory
    functions :+= scatterFunction

    // Create the gather functions for each output field
    var gatherFunctions = Map.empty[Field, GatherFunction]
    for (outputField <- originalFunction.outputFields) {

      // Create the gather function based on @Gather
      val gatherFunction = getGatherFunction(outputField)
      gatherFunction.setOriginalFunction(originalFunction)
      gatherFunction.jobNamePrefix = originalFunction.jobNamePrefix
      gatherFunction.commandDirectory = originalFunction.temp("gather-" + outputField.getName)

      val gatheredValue = originalFunction.getFieldValue(outputField).asInstanceOf[gatherFunction.GatherType]
      gatherFunction.originalOutput = gatheredValue

      tempDirectories :+= gatherFunction.commandDirectory
      cleanupFunction.originalOutputs :+= gatheredValue

      functions :+= gatherFunction

      gatherFunctions += outputField -> gatherFunction
    }

    // Create the clone functions for running the parallel jobs
    var cloneFunctions = List.empty[CommandLineFunction]
    for (i <- 1 to originalFunction.scatterCount) {
      val cloneFunction = newFunctionClone(originalFunction)
      cloneFunctions :+= cloneFunction

      val tempDir = originalFunction.temp("temp-"+i)
      cloneFunction.commandDirectory = tempDir
      tempDirectories :+= tempDir

      // Reset the input of the clone to the the temp dir and add it as an output of the scatter
      var scatterPart = CollectionUtils.updated(originalValue, resetToTempDir(tempDir))
      scatterFunction.scatterParts :+= scatterPart.asInstanceOf[scatterFunction.ScatterType]
      cloneFunction.setFieldValue(scatterField, scatterPart)

      // For each each output field, change value to the temp dir and feed it into the gatherer
      for (outputField <- originalFunction.outputFields) {
        val gatherFunction = gatherFunctions(outputField)
        val gatherPart = cloneFunction.mapField(outputField, resetToTempDir(tempDir))
        gatherFunction.gatherParts :+= gatherPart.asInstanceOf[gatherFunction.GatherType]
      }
    }
    functions = cloneFunctions ::: functions

    // Create a function to create all of the temp directories.
    // All of its inputs are the inputs of the original function.
    val initializeFunction = new CreateTempDirsFunction
    initializeFunction.jobNamePrefix = originalFunction.jobNamePrefix
    initializeFunction.commandDirectory = originalFunction.commandDirectory

    for (inputField <- originalFunction.inputFields)
      initializeFunction.originalInputs :+= originalFunction.getFieldValue(inputField)

    initializeFunction.tempDirectories = tempDirectories
    scatterFunction.tempDirectories = tempDirectories
    cleanupFunction.tempDirectories = tempDirectories

    functions +:= initializeFunction
    functions :+= cleanupFunction

    // Return all the various functions we created
    functions
  }

  private def resetToTempDir(tempDir: File): Any => Any = {
    (any: Any) => {
      any match {
        case file: File => IOUtils.reset(tempDir, file)
        case x => x
      }
    }
  }

  private def getScatterFunction(inputField: Field) =
    ReflectionUtils.getAnnotation(inputField, classOf[Scatter]).value.newInstance.asInstanceOf[ScatterFunction]

  private def getGatherFunction(outputField: Field) =
    ReflectionUtils.getAnnotation(outputField, classOf[Gather]).value.newInstance.asInstanceOf[GatherFunction]

  private def newFunctionClone(originalFunction: ScatterGatherableFunction) = {
    val cloneFunction = originalFunction.cloneFunction.asInstanceOf[ScatterGatherableFunction]
    // Make sure clone doesn't get scattered
    cloneFunction.scatterCount = 1
    cloneFunction
  }
}
