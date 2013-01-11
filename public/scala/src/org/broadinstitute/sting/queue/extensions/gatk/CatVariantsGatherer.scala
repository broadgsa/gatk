package org.broadinstitute.sting.queue.extensions.gatk

import org.broadinstitute.sting.queue.function.RetryMemoryLimit
import org.broadinstitute.sting.queue.function.scattergather.GatherFunction


/**
 * Created with IntelliJ IDEA.
 * User: ami
 * Date: 12/11/12
 * Time: 2:04 PM
 * To change this template use File | Settings | File Templates.
 */
class CatVariantsGatherer extends CatVariants with GatherFunction with RetryMemoryLimit{
  this.assumeSorted = true

  private lazy val originalGATK = this.originalFunction.asInstanceOf[CommandLineGATK]

  override def freezeFieldValues() {
    this.reference = originalGATK.reference_sequence
    this.variant = this.gatherParts.zipWithIndex map { case (input, index) => new TaggedFile(input, "input"+index) }
    this.outputFile = this.originalOutput
    super.freezeFieldValues()
  }



}
