package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.CommandLineFunction

class ShellJobManager extends JobManager[CommandLineFunction, ShellJobRunner] {
  def create(function: CommandLineFunction) = new ShellJobRunner(function)
}
