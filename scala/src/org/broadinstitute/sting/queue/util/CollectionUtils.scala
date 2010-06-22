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

  def updated(value: Any, f: (Any) => Any): Any = {
    value match {
      case seq: Seq[_] => seq.map(updated(_, f))
      case array: Array[_] => array.map(updated(_, f))
      case option: Option[_] => option.map(updated(_, f))
      case x => f(x)
    }
  }

  def foreach(value: Any, f: (Any, Any) => Unit): Unit = {
    value match {
      case seq: Seq[_] =>
        for (item <- seq) {
          f(item, seq)
          foreach(item, f)
        }
      case product: Product =>
        for (item <- product.productIterator) {
          f(item, product)
          foreach(item, f)
        }
      case array: Array[_] =>
        for (item <- array) {
          f(item, array)
          foreach(item, f)
        }
      case _ =>
    }
  }
}
