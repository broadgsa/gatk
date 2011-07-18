package org.broadinstitute.sting.queue.qscripts.examples

import org.broadinstitute.sting.queue.QScript

class HelloWorld extends QScript {
  def script() {
    add(new CommandLineFunction {
      def commandLine = "echo hello world"
    })
  }
}
