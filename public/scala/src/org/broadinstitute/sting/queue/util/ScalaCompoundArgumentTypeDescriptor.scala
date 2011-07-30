package org.broadinstitute.sting.queue.util

import collection.JavaConversions._
import org.broadinstitute.sting.queue.QException
import java.lang.Class
import org.broadinstitute.sting.commandline.{ArgumentMatches, ArgumentSource, ArgumentTypeDescriptor, ParsingEngine}
import java.lang.reflect.Type

/**
 * An ArgumentTypeDescriptor that can parse the scala collections.
 */
class ScalaCompoundArgumentTypeDescriptor extends ArgumentTypeDescriptor {

  /**
   * Checks if the class type is a scala collection.
   * @param classType Class type to check.
   * @return true if the class is a List, Set, or an Option.
   */
  def supports(classType: Class[_]) = isCompound(classType)

  /**
   * Checks if the class type is a scala collection.
   * @param source Argument source to check.
   * @return true if the source is a List, Set, or an Option.
   */
  override def isMultiValued(source: ArgumentSource) = isCompound(source.field.getType)

  /**
   * Checks if the class type is a scala collection.
   * @param classType Class type to check.
   * @return true if the class is a List, Set, or an Option.
   */
  private def isCompound(classType: Class[_]) = {
    classOf[List[_]].isAssignableFrom(classType) ||
    classOf[Set[_]].isAssignableFrom(classType) ||
    classOf[Option[_]].isAssignableFrom(classType)
  }

  /**
   * Parses the argument matches based on the class type of the argument source's field.
   * @param source Argument source that contains the field being populated.
   * @param classType Class type being parsed.
   * @param argumentMatches The argument match strings that were found for this argument source.
   * @return The parsed object.
   */
  def parse(parsingEngine: ParsingEngine, source: ArgumentSource, typeType: Type, argumentMatches: ArgumentMatches) = {
    parse(parsingEngine,source, makeRawTypeIfNecessary(typeType), argumentMatches)
  }
  
  def parse(parsingEngine: ParsingEngine, source: ArgumentSource, classType: Class[_], argumentMatches: ArgumentMatches) = {
    val componentType = ReflectionUtils.getCollectionType(source.field)
    val componentArgumentParser = parsingEngine.selectBestTypeDescriptor(componentType)

    if (classOf[List[_]].isAssignableFrom(classType)) {
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
