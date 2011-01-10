package org.broadinstitute.sting.queue.util

import org.broadinstitute.sting.BaseTest
import org.testng.annotations.Test
import java.io.File
import org.testng.Assert

class ShellJobUnitTest extends BaseTest {
  @Test
  def testEcho {
    val job = new ShellJob
    job.shellScript = writeScript("echo Hello World")
    job.run()
  }

  @Test(expectedExceptions=Array(classOf[JobExitException]))
  def testBadQuotes {
    val job = new ShellJob
    job.shellScript = writeScript("echo 'Hello World")
    job.run()
  }

  @Test
  def testGoodQuotes {
    val job = new ShellJob
    job.shellScript = writeScript("echo 'Hello World'")
    job.run()
  }

  @Test
  def testEscapeCharacters {
    var job: ShellJob = null

    job = new ShellJob
    job.shellScript = writeScript("echo #")
    job.outputFile = File.createTempFile("temp", "")
    job.run()
    Assert.assertEquals(IOUtils.readContents(job.outputFile).trim, "")

    job = new ShellJob
    job.shellScript = writeScript("""echo \#""")
    job.outputFile = File.createTempFile("temp", "")
    job.run()
    Assert.assertEquals(IOUtils.readContents(job.outputFile).trim, "#")

    job = new ShellJob
    job.shellScript = writeScript("""echo \\#""")
    job.outputFile = File.createTempFile("temp", "")
    job.run()
    Assert.assertEquals(IOUtils.readContents(job.outputFile).trim, """\#""")
  }

  @Test
  def testLongCommand {
    // This command fails on some systems with a 4096 character limit when run via the old sh -c "echo ...",
    // but works on the same systems when run via sh <script>
    val builder = new StringBuilder
    builder.append("echo ")
    for (i <- 1 to 500) {
      val s = i.toString
      builder.append("000".take(3-s.length)).append(s).append("______ ")
    }

    val job = new ShellJob
    job.shellScript = writeScript(builder.toString)
    job.run()
  }

  private val tempDir = new File(System.getProperty("java.io.tmpdir"))

  private def writeScript(contents: String) = {
    val script = IOUtils.writeTempFile(contents, "temp", "", tempDir)
    script.deleteOnExit
    script
  }
}
