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

/**
 * An importable object that provides automatic primitive to option conversion.
 */
object PrimitiveOptionConversions {

  // Conversion from Option
  implicit def byteOption2byte(x: Option[Byte]) = x.get
  implicit def shortOption2short(x: Option[Short]) = x.get
  implicit def charOption2char(x: Option[Char]) = x.get
  implicit def intOption2int(x: Option[Int]) = x.get
  implicit def longOption2long(x: Option[Long]) = x.get
  implicit def floatOption2float(x: Option[Float]) = x.get
  implicit def doubleOption2double(x: Option[Double]) = x.get
  implicit def booleanOption2boolean(x: Option[Boolean]) = x.get

  // Conversion to Option
  implicit def byte2byteOption(x: Byte) = Some(x)
  implicit def short2shortOption(x: Short) = Some(x)
  implicit def char2charOption(x: Char) = Some(x)
  implicit def int2intOption(x: Int) = Some(x)
  implicit def long2longOption(x: Long) = Some(x)
  implicit def float2floatOption(x: Float) = Some(x)
  implicit def double2doubleOption(x: Double) = Some(x)
  implicit def boolean2booleanOption(x: Boolean) = Some(x)

  // Narrowing for constants to byte, short, and float
  implicit def int2byteOption(x: Int) = Some(x.toByte)
  implicit def int2shortOption(x: Int) = Some(x.toShort)
  implicit def double2floatOption(x: Float) = Some(x.toFloat)

  // Widening
  implicit def byte2shortOption(x: Byte) = Some(x.toShort)
  implicit def byte2intOption(x: Byte) = Some(x.toInt)
  implicit def byte2longOption(x: Byte) = Some(x.toLong)
  implicit def byte2floatOption(x: Byte) = Some(x.toFloat)
  implicit def byte2doubleOption(x: Byte) = Some(x.toDouble)

  implicit def short2intOption(x: Short) = Some(x.toInt)
  implicit def short2longOption(x: Short) = Some(x.toLong)
  implicit def short2floatOption(x: Short) = Some(x.toFloat)
  implicit def short2doubleOption(x: Short) = Some(x.toDouble)

  implicit def char2intOption(x: Char) = Some(x.toInt)
  implicit def char2longOption(x: Char) = Some(x.toLong)
  implicit def char2floatOption(x: Char) = Some(x.toFloat)
  implicit def char2doubleOption(x: Char) = Some(x.toDouble)

  implicit def int2longOption(x: Int) = Some(x.toLong)
  implicit def int2floatOption(x: Int) = Some(x.toFloat)
  implicit def int2doubleOption(x: Int) = Some(x.toDouble)

  implicit def long2floatOption(x: Long) = Some(x.toFloat)
  implicit def long2doubleOption(x: Long) = Some(x.toDouble)

  implicit def float2doubleOption(x: Float) = Some(x.toDouble)

}

/**
 * A trait that exposes the above functions to all sub classes as well.
 */
trait PrimitiveOptionConversions {
  // How to we import these implicit definitions into the trait so that they are seen by objects extending a trait?
  // import PrimitiveOptionConversion._ inside of a trait does not seem to work?
  // Declaring them in a trait like this does work but does not seem scala-ish.

  implicit def byteOption2byte(x: Option[Byte]) = PrimitiveOptionConversions.byteOption2byte(x)
  implicit def shortOption2short(x: Option[Short]) = PrimitiveOptionConversions.shortOption2short(x)
  implicit def charOption2char(x: Option[Char]) = PrimitiveOptionConversions.charOption2char(x)
  implicit def intOption2int(x: Option[Int]) = PrimitiveOptionConversions.intOption2int(x)
  implicit def longOption2long(x: Option[Long]) = PrimitiveOptionConversions.longOption2long(x)
  implicit def floatOption2float(x: Option[Float]) = PrimitiveOptionConversions.floatOption2float(x)
  implicit def doubleOption2double(x: Option[Double]) = PrimitiveOptionConversions.doubleOption2double(x)
  implicit def booleanOption2boolean(x: Option[Boolean]) = PrimitiveOptionConversions.booleanOption2boolean(x)

  implicit def byte2byteOption(x: Byte) = PrimitiveOptionConversions.byte2byteOption(x)
  implicit def short2shortOption(x: Short) = PrimitiveOptionConversions.short2shortOption(x)
  implicit def char2charOption(x: Char) = PrimitiveOptionConversions.char2charOption(x)
  implicit def int2intOption(x: Int) = PrimitiveOptionConversions.int2intOption(x)
  implicit def long2longOption(x: Long) = PrimitiveOptionConversions.long2longOption(x)
  implicit def float2floatOption(x: Float) = PrimitiveOptionConversions.float2floatOption(x)
  implicit def double2doubleOption(x: Double) = PrimitiveOptionConversions.double2doubleOption(x)
  implicit def boolean2booleanOption(x: Boolean) = PrimitiveOptionConversions.boolean2booleanOption(x)

  implicit def int2byteOption(x: Int) = PrimitiveOptionConversions.int2byteOption(x)
  implicit def int2shortOption(x: Int) = PrimitiveOptionConversions.int2shortOption(x)
  implicit def double2floatOption(x: Float) = PrimitiveOptionConversions.double2floatOption(x)

  implicit def byte2shortOption(x: Byte) = PrimitiveOptionConversions.byte2shortOption(x)
  implicit def byte2intOption(x: Byte) = PrimitiveOptionConversions.byte2intOption(x)
  implicit def byte2longOption(x: Byte) = PrimitiveOptionConversions.byte2longOption(x)
  implicit def byte2floatOption(x: Byte) = PrimitiveOptionConversions.byte2floatOption(x)
  implicit def byte2doubleOption(x: Byte) = PrimitiveOptionConversions.byte2doubleOption(x)

  implicit def short2intOption(x: Short) = PrimitiveOptionConversions.short2intOption(x)
  implicit def short2longOption(x: Short) = PrimitiveOptionConversions.short2longOption(x)
  implicit def short2floatOption(x: Short) = PrimitiveOptionConversions.short2floatOption(x)
  implicit def short2doubleOption(x: Short) = PrimitiveOptionConversions.short2doubleOption(x)

  implicit def char2intOption(x: Char) = PrimitiveOptionConversions.char2intOption(x)
  implicit def char2longOption(x: Char) = PrimitiveOptionConversions.char2longOption(x)
  implicit def char2floatOption(x: Char) = PrimitiveOptionConversions.char2floatOption(x)
  implicit def char2doubleOption(x: Char) = PrimitiveOptionConversions.char2doubleOption(x)

  implicit def int2longOption(x: Int) = PrimitiveOptionConversions.int2longOption(x)
  implicit def int2floatOption(x: Int) = PrimitiveOptionConversions.int2floatOption(x)
  implicit def int2doubleOption(x: Int) = PrimitiveOptionConversions.int2doubleOption(x)

  implicit def long2floatOption(x: Long) = PrimitiveOptionConversions.long2floatOption(x)
  implicit def long2doubleOption(x: Long) = PrimitiveOptionConversions.long2doubleOption(x)

  implicit def float2doubleOption(x: Float) = PrimitiveOptionConversions.float2doubleOption(x)
}
