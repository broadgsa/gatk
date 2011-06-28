/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.queue.pipeline.examples

import org.testng.annotations.Test
import org.broadinstitute.sting.queue.pipeline.{PipelineTest, PipelineTestSpec}

class HelloWorldPipelineTest {
  @Test
  def testHelloWorld {
    val spec = new PipelineTestSpec
    spec.name = "HelloWorld"
    spec.args = "-S scala/qscript/examples/HelloWorld.scala"
    PipelineTest.executeTest(spec)
  }

  @Test
  def testHelloWorldWithPrefix {
    val spec = new PipelineTestSpec
    spec.name = "HelloWorldWithPrefix"
    spec.args = "-S scala/qscript/examples/HelloWorld.scala -jobPrefix HelloWorld"
    PipelineTest.executeTest(spec)
  }

  @Test(enabled=false)
  def testHelloWorldWithPriority {
    val spec = new PipelineTestSpec
    spec.name = "HelloWorldWithPriority"
    spec.args = "-S scala/qscript/examples/HelloWorld.scala -jobPriority 100"
    PipelineTest.executeTest(spec)
  }
}
