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

import org.broadinstitute.sting.commandline.ArgumentSource
import java.io.File
import org.broadinstitute.sting.queue.function.{QFunction, CommandLineFunction}

/**
 * Shadow clones another command line function.
 */
class CloneFunction extends CommandLineFunction {
  var originalFunction: ScatterGatherableFunction = _
  var cloneIndex: Int = _

  private var overriddenFields = Map.empty[ArgumentSource, Any]
  private var withScatterPartCount = 0

  private def withScatterPart[A](f: () => A): A = {
    var originalValues = Map.empty[ArgumentSource, Any]
    withScatterPartCount += 1
    if (withScatterPartCount == 1) {
      overriddenFields.foreach{
        case (field, overrideValue) => {
          originalValues += field -> originalFunction.getFieldValue(field)
          originalFunction.setFieldValue(field, overrideValue)
        }
      }
    }
    try {
      f()
    } finally {
      if (withScatterPartCount == 1) {
        originalValues.foreach{
          case (name, value) =>
            originalFunction.setFieldValue(name, value)
        }
      }
      withScatterPartCount -= 1
    }
  }

  override def commandOutputs = withScatterPart(() => originalFunction.commandOutputs)
  override def dotString = withScatterPart(() => originalFunction.dotString)
  override def description = withScatterPart(() => originalFunction.description)
  override protected def functionFieldClass = originalFunction.getClass

  def commandLine = withScatterPart(() => originalFunction.commandLine)

  def getFieldValue(field: String): AnyRef = {
    val source = QFunction.findField(originalFunction.getClass, field)
    getFieldValue(source)
  }

  override def getFieldValue(source: ArgumentSource): AnyRef = {
    source.field.getName match {
      case "jobOutputFile" => jobOutputFile
      case "jobErrorFile" => jobErrorFile
      case _ => overriddenFields.get(source) match {
        case Some(value) => value.asInstanceOf[AnyRef]
        case None => {
          val value = originalFunction.getFieldValue(source)
          overriddenFields += source -> value
          value
        }
      }
    }
  }

  def setFieldValue(field: String, value: Any): Unit = {
    val source = QFunction.findField(originalFunction.getClass, field)
    setFieldValue(source, value)
  }

  override def setFieldValue(source: ArgumentSource, value: Any): Unit = {
    source.field.getName match {
      case "jobOutputFile" => jobOutputFile = value.asInstanceOf[File]
      case "jobErrorFile" => jobErrorFile = value.asInstanceOf[File]
      case _ => overriddenFields += source -> value
    }
  }
}
