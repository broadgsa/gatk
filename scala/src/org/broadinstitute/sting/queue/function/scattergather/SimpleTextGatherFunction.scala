package org.broadinstitute.sting.queue.function.scattergather

import org.broadinstitute.sting.queue.function.InProcessFunction
import org.broadinstitute.sting.queue.QException
import java.io.PrintWriter
import org.apache.commons.io.{LineIterator, IOUtils, FileUtils}

/**
 * Merges a text file.
 */
class SimpleTextGatherFunction extends GatherFunction with InProcessFunction {
  def run() = {
    waitForGatherParts
    if (gatherParts.size < 1) {
      throw new QException("No files to gather to output: " + originalOutput)
    } else if (gatherParts.size == 1) {
      FileUtils.copyFile(gatherParts(0), originalOutput)
    } else {
      val writer = new PrintWriter(originalOutput)
      try {
        var startLine = 0

        val readerA = FileUtils.lineIterator(gatherParts(0))
        val readerB = FileUtils.lineIterator(gatherParts(1))
        try {
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
        } finally {
          LineIterator.closeQuietly(readerA)
          LineIterator.closeQuietly(readerB)
        }

        for (file <- gatherParts) {
          val reader = FileUtils.lineIterator(file)
          try {
            var lineNum = 0
            while (reader.hasNext && lineNum < startLine) {
              reader.nextLine
              lineNum += 1
            }
            while (reader.hasNext)
              writer.println(reader.nextLine)
          } finally {
            LineIterator.closeQuietly(reader)
          }
        }
      } finally {
        IOUtils.closeQuietly(writer)
      }
    }
  }
}
