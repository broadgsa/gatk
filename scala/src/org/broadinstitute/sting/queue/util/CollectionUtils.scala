package org.broadinstitute.sting.queue.util

/**
 * Utilities that try to deeply apply operations to collections
 */
object CollectionUtils {

  def test(value: Any, f: Any => Boolean): Boolean = {
    var result = f(value)
    foreach(value, (item, collection) => {
      result |= f(item)
    })
    result
  }

  def updated(value: Any, f: Any => Any): Any = {
    value match {
      case traversable: Traversable[_] => traversable.map(updated(_, f))
      case option: Option[_] => option.map(updated(_, f))
      case x => f(x)
    }
  }

  def foreach(value: Any, f: (Any, Any) => Unit): Unit = {
    value match {
      case traversable: Traversable[_] =>
        for (item <- traversable) {
          f(item, traversable)
          foreach(item, f)
        }
      case option: Option[_] =>
        for (item <- option) {
          f(item, option)
          foreach(item, f)
        }
      case _ =>
    }
  }

  // Because scala allows but throws NPE when trying to hash a collection with a null in it.
  // http://thread.gmane.org/gmane.comp.lang.scala.internals/3267
  // https://lampsvn.epfl.ch/trac/scala/ticket/2935
  def removeNullOrEmpty[T](value: T): T = filterNotNullOrNotEmpty(value)

  private def filterNotNullOrNotEmpty[T](value: T): T = {
    val newValue = value match {
      case traversable: Traversable[_] => traversable.map(filterNotNullOrNotEmpty(_)).filter(isNotNullOrNotEmpty(_)).asInstanceOf[T]
      case option: Option[_] => option.map(filterNotNullOrNotEmpty(_)).filter(isNotNullOrNotEmpty(_)).asInstanceOf[T]
      case x => x
    }
    newValue
  }

  private def isNotNullOrNotEmpty(value: Any): Boolean = {
    val result = value match {
      case traversable: Traversable[_] => !filterNotNullOrNotEmpty(traversable).isEmpty
      case option: Option[_] => !filterNotNullOrNotEmpty(option).isEmpty
      case x => x != null
    }
    result
  }
}
