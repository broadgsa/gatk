package org.broadinstitute.sting.queue.function.scattergather

import org.broadinstitute.sting.commandline.Argument
import org.broadinstitute.sting.queue.function.InProcessFunction
import org.apache.commons.io.FileUtils
import org.broadinstitute.sting.queue.QException
import collection.JavaConversions._
import java.io.PrintWriter

/**
 * Merges a text file.
 * The script can be changed by setting rmdirScript.
 * By default uses mergeText.sh in Sting/shell.
 * The format of the call is <mergeTextScript> <file_output> <file_1> [.. <file_n>]
 */
class SimpleTextGatherFunction extends GatherFunction with InProcessFunction {
  @Argument(doc="merge text script")
  var mergeTextScript = "mergeText.sh"


  def run() = {
    if (gatherParts.size < 1) {
      throw new QException("No files to gather to output: " + originalOutput)
    } else if (gatherParts.size == 1) {
      FileUtils.copyFile(gatherParts(0), originalOutput)
    } else {
      val writer = new PrintWriter(originalOutput)
      var startLine = 0

      val readerA = FileUtils.lineIterator(gatherParts(0))
      val readerB = FileUtils.lineIterator(gatherParts(1))
      var headersMatch = true
      while (headersMatch) {
        if (readerA.hasNext && readerB.hasNext) {
          val headerA = readerA.nextLine
          val headerB = readerB.nextLine
          headersMatch = headerA == headerB
          if (headersMatch) {
            startLine += 1
            writer.println(headerA)
          }
        } else {
          headersMatch = false
        }
      }
      readerA.close
      readerB.close

      for (file <- gatherParts) {
        val reader = FileUtils.lineIterator(file)
        var lineNum = 0
        while (reader.hasNext && lineNum < startLine) {
          reader.nextLine
          lineNum += 1
        }
        while (reader.hasNext)
          writer.println(reader.nextLine)
      }
      writer.close
    }
  }
}
