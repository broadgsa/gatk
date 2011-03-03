package org.broadinstitute.sting.queue.extensions.gatk

import org.broadinstitute.sting.queue.function.scattergather.{ScatterGatherableFunction, GatherFunction}

/**
 * Merges a vcf text file.
 */
class VcfGatherFunction extends CombineVariants with GatherFunction {

  private var originalGATK: CommandLineGATK = _

  override def setScatterGatherable(originalFunction: ScatterGatherableFunction) {
    this.originalGATK = originalFunction.asInstanceOf[CommandLineGATK]
  }

  override def freezeFieldValues = {
    this.memoryLimit = Some(1)

    this.jarFile = this.originalGATK.jarFile
    this.reference_sequence = this.originalGATK.reference_sequence
    this.intervals = this.originalGATK.intervals
    this.intervalsString = this.originalGATK.intervalsString

    this.rodBind = this.gatherParts.zipWithIndex map { case (input, index) => new RodBind("input"+index, "VCF", input) }
    this.rod_priority_list = (0 until this.gatherParts.size).map("input"+_).mkString(",")
    this.out = this.originalOutput
    this.assumeIdenticalSamples = true

    super.freezeFieldValues
  }
}
