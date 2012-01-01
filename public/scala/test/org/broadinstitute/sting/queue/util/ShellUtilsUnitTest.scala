package org.broadinstitute.sting.queue.util

import org.testng.annotations.Test
import org.testng.Assert
import java.io.{InputStreamReader, BufferedReader}

class ShellUtilsUnitTest {

  @Test
  def testEscapeShellArgumentOneCharSequences() {
    // Test all ASCII characters except \0, \n, and \r, which we do not support escaping
    for ( asciiCode <- 1 to 127 if asciiCode != 10 && asciiCode != 13  ) {
      val originalString: String = "%c".format(asciiCode.toChar)
      val quotedString: String = ShellUtils.escapeShellArgument(originalString)

      val child : Process = new ProcessBuilder("/bin/sh", "-c", "printf \"%s\" " + quotedString).start()
      val childReader : BufferedReader = new BufferedReader(new InputStreamReader(child.getInputStream))
      val childOutputBuffer : StringBuilder = new StringBuilder

      val childReaderThread : Thread = new Thread(new Runnable() {
        def run() {
          var line : String = childReader.readLine()

          while ( line != null ) {
            childOutputBuffer.append(line)
            line = childReader.readLine()
          }
        }
      })
      childReaderThread.start()

      val childReturnValue = child.waitFor()
      childReaderThread.join()

      childReader.close()
      val childOutput = childOutputBuffer.toString()

      if ( childReturnValue != 0 ) {
        Assert.fail("With character ASCII %d, sh child process returned: %d".format(asciiCode, childReturnValue))
      }
      else if ( ! originalString.equals(childOutput) ) {
        Assert.fail("With character ASCII %d, sh child process output \"%s\" instead of the expected \"%s\"".format(
                    asciiCode, childOutput, originalString))
      }
    }
  }

  @Test(expectedExceptions = Array(classOf[IllegalArgumentException]))
  def testEscapeShellArgumentNullString() {
    ShellUtils.escapeShellArgument(null)
  }

  @Test
  def testEscapeShellArgumentEmptyString() {
    Assert.assertEquals(ShellUtils.escapeShellArgument(""), "''")
  }
}