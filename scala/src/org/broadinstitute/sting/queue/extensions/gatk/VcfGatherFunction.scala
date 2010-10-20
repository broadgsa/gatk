package org.broadinstitute.sting.queue.extensions.gatk

import org.broadinstitute.sting.queue.function.InProcessFunction
import org.broadinstitute.sting.queue.QException
import org.broadinstitute.sting.queue.function.scattergather.GatherFunction
import java.io.{FileReader, PrintWriter}
import org.apache.commons.io.{LineIterator, IOUtils, FileUtils}

/**
 * Merges a vcf text file.
 */
class VcfGatherFunction extends GatherFunction with InProcessFunction {
  def run() = {
    if (gatherParts.size < 1) {
      throw new QException("No files to gather to output: " + originalOutput)
    } else {
      val writer = new PrintWriter(originalOutput)
      try {
        var reader = new FileReader(gatherParts(0))
        try {
          IOUtils.copy(reader, writer)
        } finally {
          IOUtils.closeQuietly(reader)
        }

        for (file <- gatherParts.tail) {
          var inHeaders = true
          val itor = FileUtils.lineIterator(file)
          try {
            while (itor.hasNext) {
              val nextLine = itor.nextLine
              if (inHeaders && nextLine(0) != '#')
                inHeaders = false
              if (!inHeaders)
                writer.println(nextLine)
            }
          } finally {
            LineIterator.closeQuietly(itor)
          }
        }
      } finally {
        IOUtils.closeQuietly(writer)
      }
    }
  }
}