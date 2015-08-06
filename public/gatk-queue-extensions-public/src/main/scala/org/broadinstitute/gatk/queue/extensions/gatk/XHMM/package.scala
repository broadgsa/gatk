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

package org.broadinstitute.gatk.queue.extensions.gatk

import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk.DoC._
import org.broadinstitute.gatk.utils.commandline._
import java.io.{File, PrintStream, PrintWriter}
import collection.JavaConversions._
import org.broadinstitute.gatk.queue.function.scattergather.{CloneFunction, ScatterFunction, GatherFunction, ScatterGatherableFunction}
import org.broadinstitute.gatk.queue.function.{CommandLineFunction, InProcessFunction}
import org.broadinstitute.gatk.utils.io.IOUtils

// Minimal refactor from a package object to a file full of classes/objects
// due to ongoing bugs with inner classes/objects in package objects:
//   https://issues.scala-lang.org/browse/SI-4344
//   https://issues.scala-lang.org/browse/SI-5954

    abstract class BaseGenotypeCNVs(inputParam: File, xcnv: File, origRDParam: File, outName: String, xhmmParamsArg: File, referenceFile: File, genotypeCommandLineParams: String, xhmmExec: File, groups: List[Group]) extends SamplesScatterable(xhmmExec, groups) {
      @Input(doc = "")
      val input = inputParam

      @Input(doc = "")
      val xhmmParams = xhmmParamsArg

      @Input(doc = "")
      val origRD = origRDParam

      @Input(doc = "")
      val inXcnv = xcnv

      @Output
      @Gather(classOf[MergeVCFsGatherFunction])
      val vcf: File = new File(outName)

      override def commandLine =
        xhmmExec + " --genotype" +
          " -p " + xhmmParams +
          " -r " + input +
          " -g " + inXcnv +
          " -F " + referenceFile +
          " -R " + origRD +
          " -v " +  vcf +
          " " + genotypeCommandLineParams +
          " " + addCommand
    }

    class GenotypeCNVs(inputParam: File, xcnv: File, origRDParam: File, genotypeOutputBase: File, xhmmParamsArg: File, referenceFile: File, genotypeCommandLineParams: String, xhmmExec: File, groups: List[Group]) extends BaseGenotypeCNVs(inputParam, xcnv, origRDParam, genotypeOutputBase.getPath + ".vcf", xhmmParamsArg, referenceFile, genotypeCommandLineParams, xhmmExec, groups) {
      override def description = "Genotypes CNV regions in all samples: " + commandLine
    }


abstract class SamplesScatterable(val xhmmExec: File, val groups: List[Group]) extends ScatterGatherableFunction with CommandLineFunction {
  this.scatterCount = groups.size
  this.scatterClass = classOf[SamplesScatterFunction]

  @Input(doc = "", required=false)
  var keepSampleIDs: Option[String] = None

  def addCommand = if (keepSampleIDs.isDefined) ("--keepSampleIDs " + keepSampleIDs.get) else ""
}

class SamplesScatterFunction extends ScatterFunction with InProcessFunction {
  protected var groups: List[Group] = _
  override def scatterCount = groups.size

  @Output(doc="Scatter function outputs")
  var scatterSamples: Seq[File] = Nil

  override def init() {
    this.groups = this.originalFunction.asInstanceOf[SamplesScatterable].groups
  }

  override def bindCloneInputs(cloneFunction: CloneFunction, index: Int) {
    val scatterPart = IOUtils.absolute(cloneFunction.commandDirectory, "keepSampleIDs.txt")
    cloneFunction.setFieldValue("keepSampleIDs", Some(scatterPart))
    this.scatterSamples :+= scatterPart
  }

  override def run() {
    if (groups.size != this.scatterSamples.size)
      throw new Exception("Internal inconsistency error in scattering jobs")

    (groups, this.scatterSamples).zipped foreach {
      (group, sampsFile) => {
        val sampsWriter = new PrintWriter(new PrintStream(sampsFile))

        for (samp <- group.samples) {
          try {
            sampsWriter.printf("%s%n", samp)
          }
          catch {
            case e: Exception => throw e
          }
        }
        sampsWriter.close
      }
    }
  }
}

trait MergeVCFs extends CommandLineFunction {
  var xhmmExec: File = _

  @Input(doc = "")
  var inputVCFs: List[File] = Nil

  @Output
  var mergedVCF: File = null

  override def commandLine =
    xhmmExec + " --mergeVCFs" +
      inputVCFs.map(input => " --mergeVCF " + input).reduceLeft(_ + "" + _) +
      " -v " + mergedVCF

  override def description = "Combines VCF outputs for multiple samples (at same loci): " + commandLine
}

class MergeVCFsGatherFunction extends MergeVCFs with GatherFunction {
  override def freezeFieldValues() {
    super.freezeFieldValues()

    this.xhmmExec = originalFunction.asInstanceOf[SamplesScatterable].xhmmExec

    this.inputVCFs = this.gatherParts.toList
    this.mergedVCF = this.originalOutput
  }
}

class DummyGatherFunction extends InProcessFunction with GatherFunction {
  override def run() {}
}
