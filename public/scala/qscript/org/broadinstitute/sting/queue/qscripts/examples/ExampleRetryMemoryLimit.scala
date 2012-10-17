import org.broadinstitute.sting.queue.function.RetryMemoryLimit
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

class ExampleRetryMemoryLimit extends QScript {
  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _

  @Input(doc="Bam file to genotype.", shortName="I")
  var bamFile: File = _

  def script() {
    for (scatterCount <- 1 to 2) {
      val ug = new UnifiedGenotyper with RetryMemoryLimit
      // First run with 1m
      ug.memoryLimit = .001
      // On retry run with 1g
      ug.retryMemoryFunction = (d => d * 1000)
      ug.reference_sequence = referenceFile
      ug.input_file = Seq(bamFile)
      ug.out = swapExt(bamFile, ".bam", ".scattered_%d.vcf".format(scatterCount))
      ug.scatterCount = scatterCount
      add(ug)
    }
  }
}
