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

import collection.JavaConversions._
import org.broadinstitute.gatk.queue.QException
import java.lang.Class
import org.broadinstitute.gatk.utils.commandline.{ArgumentMatches, ArgumentSource, ArgumentTypeDescriptor, ParsingEngine}
import org.broadinstitute.gatk.utils.exceptions.UserException
import java.lang.reflect.Type

/**
 * An ArgumentTypeDescriptor that can parse the scala collections.
 */
class ScalaCompoundArgumentTypeDescriptor extends ArgumentTypeDescriptor {

  /**
   * Checks if the class type is a scala collection.
   * @param classType Class type to check.
   * @return true if the class is a Seq, Set, or an Option.
   */
  def supports(classType: Class[_]) = isCompound(classType)

  /**
   * Checks if the class type is a scala collection.
   * @param source Argument source to check.
   * @return true if the source is a Seq, Set, or an Option.
   */
  override def isMultiValued(source: ArgumentSource) = isCompound(source.field.getType)

  /**
   * Checks if the class type is a scala collection.
   * @param classType Class type to check.
   * @return true if the class is a Seq, Set, or an Option.
   */
  private def isCompound(classType: Class[_]) = {
    classOf[Seq[_]].isAssignableFrom(classType) ||
    classOf[List[_]].isAssignableFrom(classType) || // see comment below re: List vs. Seq
    classOf[Set[_]].isAssignableFrom(classType) ||
    classOf[Option[_]].isAssignableFrom(classType)
  }

  /**
   * Parses the argument matches based on the class type of the argument source's field.
   * @param parsingEngine Parsing engine.
   * @param source Argument source that contains the field being populated.
   * @param typeType Type of the argument source's field.
   * @param argumentMatches The argument match strings that were found for this argument source.
   * @return The parsed object.
   */
  def parse(parsingEngine: ParsingEngine, source: ArgumentSource, typeType: Type, argumentMatches: ArgumentMatches) = {
    parse(parsingEngine,source, ArgumentTypeDescriptor.makeRawTypeIfNecessary(typeType), argumentMatches)
  }
  
  def parse(parsingEngine: ParsingEngine, source: ArgumentSource, classType: Class[_], argumentMatches: ArgumentMatches) = {
    val componentType = ReflectionUtils.getCollectionType(source.field)
    if (componentType == classOf[java.lang.Object])
      throw new UserException.CannotExecuteQScript("Please also include a @ClassType(classOf[<primitive type>]) annotation on field: " + source.field + ". Example: @ClassType(classOf[Double]). The scala generic type for the field was subjected to java/scala type erasure and is not available via reflection.")
    val componentArgumentParser = parsingEngine.selectBestTypeDescriptor(componentType)

    if (classOf[Seq[_]].isAssignableFrom(classType)) {
      var seq = Seq.empty[Any]
      for (argumentMatch <- argumentMatches)
        for (value <- argumentMatch)
          seq :+= componentArgumentParser.parse(parsingEngine, source, componentType, new ArgumentMatches(value))
      seq
    } else if (classOf[List[_]].isAssignableFrom(classType)) {
      // QScripts should be using the interface Seq instead of the class List.
      // Leaving this here for now for legacy support until the effects of switching have been tested for a while. -ks
      var list = List.empty[Any]
      for (argumentMatch <- argumentMatches)
        for (value <- argumentMatch)
          list :+= componentArgumentParser.parse(parsingEngine, source, componentType, new ArgumentMatches(value))
      list
    } else if (classOf[Set[_]].isAssignableFrom(classType)) {
      var set = Set.empty[Any]
      for (argumentMatch <- argumentMatches)
        for (value <- argumentMatch)
            set += componentArgumentParser.parse(parsingEngine, source, componentType, new ArgumentMatches(value))
      set
    } else if (classOf[Option[_]].isAssignableFrom(classType)) {
      if (argumentMatches.size > 1)
        throw new QException("Unable to set Option to multiple values: " + argumentMatches.mkString(" "))
      else if (argumentMatches.size == 1)
        Some(componentArgumentParser.parse(parsingEngine, source, componentType, argumentMatches))
      else
        None
    } else
      throw new QException("Unsupported compound argument type: " + classType)
  }
}
