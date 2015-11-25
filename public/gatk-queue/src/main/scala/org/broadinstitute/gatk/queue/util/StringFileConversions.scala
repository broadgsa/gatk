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

package org.broadinstitute.gatk.queue.util

import java.io.{Serializable, File}
import scala.language.implicitConversions

/**
 * Converts String to/from File
 * The method signatures are based on the compiler errors reported by StringFileConversionsUnitTest.
 * The implementations are based on the runtime errors reported by StringFileConversionsUnitTest.
 */
object StringFileConversions {
  implicit def stringAsFile(x: String): File = {
    if (x == null) null else new File(x)
  }

  implicit def fileAsString(x: File): String = {
    if (x == null) null else x.getPath
  }

  // Possible to get the invariants, covariants, contravariants, upper type bounds, lower type bounds
  // and mixins all correct so this doesn't have to be duplicated with concrete implementations?
  // http://programming-scala.labs.oreilly.com/ch12.html is your friend.

  implicit def stringsAsFiles(x: Seq[Comparable[_ >: String with File <: Comparable[_ >: String with File <: Serializable] with Serializable] with Serializable]): Seq[File] = {
    x.map(_ match {
      case string: String => stringAsFile(string)
      case file: File => file
      case null => null
    })
  }

  implicit def filesAsStrings(x: Seq[Comparable[_ >: File with String <: Comparable[_ >: File with String <: Serializable] with Serializable] with Serializable]): Seq[String] = {
    x.map(_ match {
      case file: File => fileAsString(file)
      case string: String => string
      case null => null
    })
  }

  implicit def stringsAsFilesList(x: List[Comparable[_ >: String with File <: Comparable[_ >: String with File <: Serializable] with Serializable] with Serializable]): List[File] = {
    x.map(_ match {
      case string: String => stringAsFile(string)
      case file: File => file
      case null => null
    })
  }

  implicit def filesAsStringsList(x: List[Comparable[_ >: File with String <: Comparable[_ >: File with String <: Serializable] with Serializable] with Serializable]): List[String] = {
    x.map(_ match {
      case file: File => fileAsString(file)
      case string: String => string
      case null => null
    })
  }

}

/**
 * Converts String to/from File
 * The method signatures are based on the compiler errors reported by StringFileConversionsUnitTest.
 * The implementations are based on the runtime errors reported by StringFileConversionsUnitTest.
 */
trait StringFileConversions {
  implicit def stringAsFile(x: String): File = {
    StringFileConversions.stringAsFile(x)
  }

  implicit def fileAsString(x: File): String = {
    StringFileConversions.fileAsString(x)
  }

  implicit def stringsAsFiles(x: Seq[Comparable[_ >: String with File <: Comparable[_ >: String with File <: Serializable] with Serializable] with Serializable]): Seq[File] = {
    StringFileConversions.stringsAsFiles(x)
  }

  implicit def filesAsStrings(x: Seq[Comparable[_ >: File with String <: Comparable[_ >: File with String <: Serializable] with Serializable] with Serializable]): Seq[String] = {
    StringFileConversions.filesAsStrings(x)
  }

  implicit def stringsAsFilesList(x: List[Comparable[_ >: String with File <: Comparable[_ >: String with File <: Serializable] with Serializable] with Serializable]): List[File] = {
    StringFileConversions.stringsAsFilesList(x)
  }

  implicit def filesAsStringsList(x: List[Comparable[_ >: File with String <: Comparable[_ >: File with String <: Serializable] with Serializable] with Serializable]): List[String] = {
    StringFileConversions.filesAsStringsList(x)
  }

}
