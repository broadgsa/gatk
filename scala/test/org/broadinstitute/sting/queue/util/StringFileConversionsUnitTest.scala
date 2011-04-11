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

package org.broadinstitute.sting.queue.util

import org.testng.annotations.Test
import java.io.File
import org.testng.Assert
import StringFileConversions._

class StringFileConversionsUnitTest {
  @Test
  def testStringToFile() {
    val file = new File("foo")
    val string: String = file
    Assert.assertEquals(string, "foo")
  }

  @Test
  def testFileToString() {
    val string = "foo"
    val file: File = string
    Assert.assertEquals(file, new File("foo"))
  }

  @Test
  def testStringToFileList() {
    var files = List(new File("foo"))
    files :+= "bar"
    Assert.assertEquals(files, List(new File("foo"), new File("bar")))
  }

  @Test
  def testFileToStringList() {
    var strings = List("foo")
    strings :+= new File("bar")
    Assert.assertEquals(strings, List("foo", "bar"))
  }

  @Test
  def testStringToFileSet() {
    var files = Set(new File("foo"))
    files += "bar"
    Assert.assertEquals(files, Set(new File("foo"), new File("bar")))
  }

  @Test
  def testFileToStringSet() {
    var strings = Set("foo")
    strings += new File("bar")
    Assert.assertEquals(strings, Set("foo", "bar"))
  }

  @Test
  def testStringListToFileList() {
    var files = List(new File("foo"))
    files ++= List("bar")
    Assert.assertEquals(files, List(new File("foo"), new File("bar")))
  }

  @Test
  def testFileListToStringList() {
    var strings = List("foo")
    strings ++= List(new File("bar"))
    Assert.assertEquals(strings, List("foo", "bar"))
  }

  @Test
  def testStringSetToFileSet() {
    var files = Set(new File("foo"))
    files ++= Set("bar")
    Assert.assertEquals(files, Set(new File("foo"), new File("bar")))
  }

  @Test
  def testFileSetToStringSet() {
    var strings = Set("foo")
    strings ++= Set(new File("bar"))
    Assert.assertEquals(strings, Set("foo", "bar"))
  }
}
