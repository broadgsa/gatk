package org.broadinstitute.sting.queue.library.clf.vcf

import java.io.File
import collection.JavaConversions._
import org.broadinstitute.sting.commandline._
import org.broadinstitute.sting.queue.function.CommandLineFunction
import org.broadinstitute.sting.utils.text.XReadLines

class VCFExtractSamples(inVCF: File, outVCF: File, samples: List[String]) extends CommandLineFunction {
  @Input(doc="input VCF from which to extract samples") var inputVCF : File = inVCF
  @Output(doc="output VCF to write extracted samples to") var outputVCF : File = outVCF
  @Argument(doc="List of samples to extract from the VCF") var sampleList : List[String] = samples

  var sampleGrep : String = _

  def this(in: File, out: File, samples: File) = this(in,out, (new XReadLines(samples)).readLines.toList)

  override def freezeFieldValues = {
    this.logger.warn("Note: Using VCFExtractSamples invalidates AC/AF/AN annotations. This is an explicit warning.")
    sampleGrep = "'" + sampleList.reduceLeft(_ + "|" + _) + "'"
    super.freezeFieldValues
  }

  def commandLine = {

    var first : String = "head -n 500 %s | grep \\\\#\\\\# > %s".format(inputVCF.getAbsolutePath,outputVCF.getAbsolutePath)
    var second : String = "head -n 500 %s | grep \\\\#CHR | tr '\\t' '\\n' | awk '{print ++count\"\\t\"$1}' ".format(inputVCF.getAbsolutePath)
    second += "| egrep %s | awk '{print $1}' | tr '\\n' ',' | xargs -i cut -f1-9,\\{\\} %s | grep -v \\\\#\\\\# >> %s".format(sampleGrep,inputVCF.getAbsolutePath,outputVCF.getAbsolutePath)

    first+" ; "+second
  }

}