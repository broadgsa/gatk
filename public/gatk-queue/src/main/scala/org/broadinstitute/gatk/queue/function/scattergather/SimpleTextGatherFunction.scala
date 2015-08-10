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

package org.broadinstitute.gatk.queue.function.scattergather

import org.broadinstitute.gatk.queue.function.InProcessFunction
import org.broadinstitute.gatk.queue.QException
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
