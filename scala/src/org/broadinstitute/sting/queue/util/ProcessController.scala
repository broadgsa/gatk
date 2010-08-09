package org.broadinstitute.sting.queue.util

import java.io._
import scala.collection.mutable.{HashSet, ListMap}

/**
 * Facade to Runtime.exec() and java.lang.Process.  Handles
 * running a process to completion and returns stdout and stderr
 * as strings.  Creates separate threads for reading stdout and stderr,
 * then reuses those threads for each process most efficient use is
 * to create one of these and use it repeatedly.  Instances are not
 * thread-safe, however.
 *
 * @author originally by Michael Koehrsen ported to scala and enhanced by Khalid Shakir
 */
class ProcessController extends Logging {

  // Threads that capture stdout and stderr
  private val stdoutCapture = new OutputCapture(ProcessController.STDOUT_KEY)
  private val stderrCapture = new OutputCapture(ProcessController.STDERR_KEY)

  // Communication channels with output capture threads
  /** Holds the stdout and stderr sent to the background capture threads */
  private val toCapture = new ListMap[String, ProcessController.CapturedStreamOutput]

  /** Holds the results of the capture from the background capture threads.
   * May be the content via toCapture or an EmptyStreamOutput if the capture was interrupted. */
  private val fromCapture = new ListMap[String, ProcessController.StreamOutput]

  // Start the background threads for this controller.
  stdoutCapture.start()
  stderrCapture.start()

  /**
   * Executes                                                                                                                           a command line program with the settings and waits for it to return, processing the output on a background thread.
   * @param settings Settings to be run.
   * @return The output of the command.
   */
  def exec(settings: ProcessController.ProcessSettings): ProcessController.ProcessOutput = {
    var builder = new ProcessBuilder(settings.cmdarray:_*)
    builder.directory(settings.directory)

    if (settings.environment != null) {
      val builderEnvironment = builder.environment
      builderEnvironment.clear()
      settings.environment.foreach{case (name, value) => builderEnvironment.put(name, value)}
    }

    builder.redirectErrorStream(settings.redirectErrorStream)

    var stdout: ProcessController.StreamOutput = null
    var stderr: ProcessController.StreamOutput = null
    val process = builder.start

    ProcessController.running.add(process)
    try {
      val stdoutSettings = if (settings.stdoutSettings == null) ProcessController.EmptyStreamSettings else settings.stdoutSettings
      val stderrSettings = if (settings.stderrSettings == null) ProcessController.EmptyStreamSettings else settings.stderrSettings

      toCapture.synchronized {
        toCapture.put(ProcessController.STDOUT_KEY, new ProcessController.CapturedStreamOutput(process.getInputStream, stdoutSettings))
        toCapture.put(ProcessController.STDERR_KEY, new ProcessController.CapturedStreamOutput(process.getErrorStream, stderrSettings))
        toCapture.notifyAll()
      }

      if (settings.stdinSettings.input != null) {
        val writer = new OutputStreamWriter(process.getOutputStream)
        writer.write(settings.stdinSettings.input)
        writer.flush()
      }
      if (settings.stdinSettings.inputFile != null) {
        val reader = new FileReader(settings.stdinSettings.inputFile)
        val writer = new OutputStreamWriter(process.getOutputStream)
        val buf = new Array[Char](4096)
        var readCount = 0
        while ({readCount = reader.read(buf); readCount} >= 0)
          writer.write(buf, 0, readCount)
        writer.flush()
        reader.close()
      }
      
      try {
        process.getOutputStream.close()
        process.waitFor()
      } finally {
        while (stdout == null || stderr == null) {
          fromCapture.synchronized {
            fromCapture.remove(ProcessController.STDOUT_KEY) match {
              case Some(stream) => stdout = stream
              case None => /* ignore */
            }
            fromCapture.remove(ProcessController.STDERR_KEY) match {
              case Some(stream) => stderr = stream
              case None => /* ignore */
            }

            try {
              if (stdout == null || stderr == null)
                fromCapture.wait()
            } catch {
              case e: InterruptedException =>
                logger.error(e)
            }
          }
        }
      }
    } finally {
      ProcessController.running.remove(process)
    }

    new ProcessController.ProcessOutput(process.exitValue, stdout, stderr)
  }

  /** Ensures that the threads used to manipulate the IO for the process are cleaned up properly. */
  def close() = {
    try {
      stdoutCapture.interrupt()
      stderrCapture.interrupt()
    } catch {
      case e =>
        logger.error(e)
    }
  }

  /** calls close() */
  override def finalize = close()

  /**
   * Reads in the output of a stream on a background thread to keep the output pipe from backing up and freezing the called process.
   * @param key The stdout or stderr key for this output capture.
   */
  private class OutputCapture(private val key: String)
          extends Thread("OutputCapture-" + key + "-" + Thread.currentThread.getName) {

    setDaemon(true)

    /** Runs the capture. */
    override def run = {
      var break = false
      while (!break) {
        var processStream: ProcessController.StreamOutput = ProcessController.EmptyStreamOutput
        try {
          // Wait for a new input stream to be passed from this process controller.
          var capturedProcessStream: ProcessController.CapturedStreamOutput = null
          while (capturedProcessStream == null) {
            toCapture.synchronized {
              toCapture.remove(key) match {
                case Some(stream) => capturedProcessStream = stream
                case None => toCapture.wait()
              }
            }
          }
          // Read in the input stream
          processStream = capturedProcessStream
          capturedProcessStream.read
        } catch {
          case e: InterruptedException => {
            logger.info("OutputReader interrupted, exiting")
            break = true
          }
          case e: IOException => {
            logger.error("Error reading process output", e)
          }
        } finally {
          // Send the string back to the process controller.
          fromCapture.synchronized {
            fromCapture.put(key, processStream)
            fromCapture.notify()
          }
        }
      }
    }
  }
}

/**
 * Facade to Runtime.exec() and java.lang.Process.  Handles
 * running a process to completion and returns stdout and stderr
 * as strings.  Creates separate threads for reading stdout and stderr,
 * then reuses those threads for each process most efficient use is
 * to create one of these and use it repeatedly.  Instances are not
 * thread-safe, however.
 *
 * @author originally by Michael Koehrsen ported to scala and enhanced by Khalid Shakir
 */
object ProcessController extends Logging {

  /**
   * Settings that define how to run a process.
   * @param cmdarray Command line to run.
   * @param environment Environment settings to override System.getEnv, or null to use System.getEnv.
   * @param directory The directory to run the command in, or null to run in the current directory.
   * @param stdinSettings Settings for writing to the process stdin.
   * @param stdoutSettings Settings for capturing the process stdout.
   * @param stderrSettings Setting for capturing the process stderr.
   * @param redirectErrorStream true if stderr should be sent to stdout.
   */
  class ProcessSettings(val cmdarray: Array[String], val environment: Map[String, String], val directory: File,
                        val stdinSettings: InputStreamSettings, val stdoutSettings: OutputStreamSettings,
                        val stderrSettings: OutputStreamSettings, val redirectErrorStream: Boolean)

  /**
   * Settings that define text to write to the process stdin.
   * @param input String to write to stdin.
   * @param inputFile File to write to stdin.
   */
  class InputStreamSettings(val input: String, val inputFile: File)

  /**
   * Settings that define text to capture from a process stream.
   * @param stringSize The number of characters to capture, or -1 for unlimited.
   * @param outputFile The file to write output to, or null to skip output.
   * @param outputFileAppend true if the output file should be appended to.
   */
  class OutputStreamSettings(val stringSize: Int, val outputFile: File, val outputFileAppend: Boolean)

  /**
   * The output of a process.
   * @param exitValue The exit value.
   * @param stdout The capture of stdout as defined by the stdout OutputStreamSettings.
   * @param stderr The capture of stderr as defined by the stderr OutputStreamSettings.
   */
  class ProcessOutput(val exitValue: Int, val stdout: StreamOutput, val stderr: StreamOutput)

  /**
   * The base class of stream output.
   */
  abstract class StreamOutput {
    /**
     * Returns the content as a string.
     * @return The content as a string.
     */
    def content: String

    /**
     * Returns true if the content was truncated.
     * @return true if the content was truncated.
     */
    def contentTruncated: Boolean
  }

  private var currentCaptureId = 0
  /**
   * Returns the next output capture id.
   * @return The next output capture id.
   */
  private def NEXT_OUTPUT_CAPTURE_ID = {
    currentCaptureId += 1
    currentCaptureId
  }
  private val STDOUT_KEY = "stdout"
  private val STDERR_KEY = "stderr"

  /** Tracks running processes so that they can be killed as the JVM shuts down. */
  private val running = new HashSet[Process]()
  Runtime.getRuntime.addShutdownHook(new Thread {
    /** Kills running processes as the JVM shuts down. */
    override def run = for (process <- running.clone) {
      logger.warn("Killing: " + process)
      process.destroy
    }
  })

  /** Empty stream settings used when no output is requested. */
  private object EmptyStreamSettings extends OutputStreamSettings(0, null, false)

  /** Empty stream output when no output is captured due to an error. */
  private object EmptyStreamOutput extends StreamOutput {
    def content = ""
    def contentTruncated = false
  }

  /**
   * Stream output captured from a stream.
   * @param stream Stream to capture output.
   * @param settings Settings that define what to capture.
   */
  private class CapturedStreamOutput(val stream: InputStream, val settings: OutputStreamSettings) extends StreamOutput {
    /**
     * Returns the captured content as a string.
     * @return The captured content as a string.
     */
    def content = stringWriter.toString()

    /**
     * Returns true if the captured content was truncated.
     * @return true if the captured content was truncated.
     */
    def contentTruncated = stringTruncated

    /**
     * Drain the input stream to keep the process from backing up until it's empty.
     */
    def read() = {
      val reader = new InputStreamReader(stream)
      val buf = new Array[Char](4096)
      var readCount = 0
      while ({readCount = reader.read(buf); readCount} >= 0) {
        writeString(buf, readCount)
        writeFile(buf, readCount)
      }
      closeFile()
      stream.close()
    }

    /** The string to write capture content. */
    private lazy val stringWriter = if (settings.stringSize < 0) new StringWriter else new StringWriter(settings.stringSize)

    /** True if the content is truncated. */
    private var stringTruncated = false

    /** The number of characters left until the buffer is full. */
    private var stringRemaining = settings.stringSize

    /**
     * Writes the buffer to the stringWriter up to stringRemaining characters.
     * @param chars Character buffer to write.
     * @param len Number of characters in the buffer.
     */
    private def writeString(chars: Array[Char], len: Int) = {
      if (settings.stringSize < 0) {
        stringWriter.write(chars, 0, len)
      } else {
        if (!stringTruncated) {
          stringWriter.write(chars, 0, if (len > stringRemaining) stringRemaining else len)
          stringRemaining -= len
          if (stringRemaining < 0)
            stringTruncated = true
        }
      }
    }

    /** The file writer to capture content or null if no output file was requested. */
    private lazy val fileWriter = {
      if (settings.outputFile == null) {
        null
      } else {
        new FileWriter(settings.outputFile, settings.outputFileAppend)
      }
    }

    /**
     * Writes the buffer to the fileWriter if it is not null.
     * @param chars Character buffer to write.
     * @param len Number of characters in the buffer.
     */
    private def writeFile(chars: Array[Char], len: Int) = {
      if (fileWriter != null)
        fileWriter.write(chars, 0, len)
    }

    /** Closes the fileWriter if it is not null. */
    private def closeFile() = {
      if (fileWriter != null) {
        fileWriter.flush
        fileWriter.close
      }
    }
  }
}
