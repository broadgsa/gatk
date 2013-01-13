/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.queue.util

import org.testng.annotations.Test
import java.io.File
import org.testng.Assert
import StringFileConversions._

class StringFileConversionsUnitTest {
  @Test
  def testStringToFile() {
    var string: String = new File("foo")
    Assert.assertEquals(string, "foo")

    string = null.asInstanceOf[File]
    Assert.assertNull(string)
  }

  @Test
  def testFileToString() {
    var file: File = "foo"
    Assert.assertEquals(file, new File("foo"))

    file = null.asInstanceOf[String]
    Assert.assertNull(file)
  }

  @Test
  def testStringToFileList() {
    var files = Seq(new File("foo"))
    files :+= "bar"
    Assert.assertEquals(files, Seq(new File("foo"), new File("bar")))

    files = Seq(new File("foo"))
    files :+= null.asInstanceOf[String]
    Assert.assertEquals(files, Seq(new File("foo"), null))

    files = Seq[File](null)
    files :+= "foo"
    Assert.assertEquals(files, Seq(null, new File("foo")))

    files = Seq[File](null)
    files :+= null.asInstanceOf[String]
    Assert.assertEquals(files, Seq(null, null))
  }

  @Test
  def testFileToStringList() {
    var strings = Seq("foo")
    strings :+= new File("bar")
    Assert.assertEquals(strings, Seq("foo", "bar"))

    strings = Seq("foo")
    strings :+= null.asInstanceOf[File]
    Assert.assertEquals(strings, Seq("foo", null))

    strings = Seq[String](null)
    strings :+= new File("foo")
    Assert.assertEquals(strings, Seq(null, "foo"))

    strings = Seq[String](null)
    strings :+= null.asInstanceOf[File]
    Assert.assertEquals(strings, Seq(null, null))
  }

  @Test
  def testStringToFileSet() {
    var files = Set(new File("foo"))
    files += "bar"
    Assert.assertEquals(files, Set(new File("foo"), new File("bar")))

    files = Set(new File("foo"))
    files += null.asInstanceOf[String]
    Assert.assertEquals(files, Set(new File("foo"), null))

    files = Set[File](null)
    files += "foo"
    Assert.assertEquals(files, Set(new File("foo"), null))

    files = Set[File](null)
    files += null.asInstanceOf[String]
    Assert.assertEquals(files, Set(null))
  }

  @Test
  def testFileToStringSet() {
    var strings = Set("foo")
    strings += new File("bar")
    Assert.assertEquals(strings, Set("foo", "bar"))

    strings = Set("foo")
    strings += null.asInstanceOf[File]
    Assert.assertEquals(strings, Set("foo", null))

    strings = Set[String](null)
    strings += new File("foo")
    Assert.assertEquals(strings, Set("foo", null))

    strings = Set[String](null)
    strings += null.asInstanceOf[File]
    Assert.assertEquals(strings, Set(null))
  }

  @Test
  def testStringListToFileList() {
    var files = Seq(new File("foo"))
    files ++= Seq("bar")
    Assert.assertEquals(files, Seq(new File("foo"), new File("bar")))

    files = Seq(new File("foo"))
    files ++= Seq[String](null)
    Assert.assertEquals(files, Seq(new File("foo"), null))

    files = Seq[File](null)
    files ++= Seq("foo")
    Assert.assertEquals(files, Seq(null, new File("foo")))

    files = Seq[File](null)
    files ++= Seq[String](null)
    Assert.assertEquals(files, Seq(null, null))
  }

  @Test
  def testFileListToStringList() {
    var strings = Seq("foo")
    strings ++= Seq(new File("bar"))
    Assert.assertEquals(strings, Seq("foo", "bar"))

    strings = Seq("foo")
    strings ++= Seq[File](null)
    Assert.assertEquals(strings, Seq("foo", null))

    strings = Seq[String](null)
    strings ++= Seq(new File("foo"))
    Assert.assertEquals(strings, Seq(null, "foo"))

    strings = Seq[String](null)
    strings ++= Seq[File](null)
    Assert.assertEquals(strings, Seq(null, null))
  }

  @Test
  def testStringSetToFileSet() {
    var files = Set(new File("foo"))
    files ++= Set("bar")
    Assert.assertEquals(files, Set(new File("foo"), new File("bar")))

    files = Set(new File("foo"))
    files ++= Set[String](null)
    Assert.assertEquals(files, Set(new File("foo"), null))

    files = Set[File](null)
    files ++= Set("foo")
    Assert.assertEquals(files, Set(new File("foo"), null))

    files = Set[File](null)
    files ++= Set[String](null)
    Assert.assertEquals(files, Set(null))
  }

  @Test
  def testFileSetToStringSet() {
    var strings = Set("foo")
    strings ++= Set(new File("bar"))
    Assert.assertEquals(strings, Set("foo", "bar"))

    strings = Set("foo")
    strings ++= Set[File](null)
    Assert.assertEquals(strings, Set("foo", null))

    strings = Set[String](null)
    strings ++= Set(new File("foo"))
    Assert.assertEquals(strings, Set("foo", null))

    strings = Set[String](null)
    strings ++= Set[File](null)
    Assert.assertEquals(strings, Set(null))
  }
}
