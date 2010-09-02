package modiphy.util

trait Memo[A,B] extends Function[A,B]{
  def apply(a:A):B
}


class ArrayMemo[A,B <: AnyRef](g:A=>Int,size:Int)(f:A=>B)(implicit m:Manifest[B]) extends Memo[A,B]{
  val memo= new Array[B](size)
  var junk = List[B]()
  def apply(a:A)={
    val num = g(a)
    memo(num) match {
      case null => val ans = f(a); memo(num)=ans; ans
      case ans => ans
    }
  }
}
class PermMemo[A, B](f:A=>B) extends Memo[A,B]{
  var map=Map[A,B]()
  def apply(a:A)={
    if (map contains a){map(a)}
    else {val ans = f(a);
    map=map updated (a,ans)
    ans}
  }
}

class ConcurrentPermMemo[A, B](f:A=>B) extends Memo[A,B]{
  var map=Map[A,B]()
  def apply(a:A)=Symbol(a.hashCode.toString).synchronized{
    if (map contains a){map(a)}
    else {val ans = f(a);
    map.synchronized{map=map updated (a,ans)}
    ans}
  }
}
class LRUMemo[A,B](size:Int,f:A=>B) extends Memo[A,B]{
  val map=new CacheMap[A,B](size)
  def apply(a:A)={
    if (map contains a){map(a)}
    else {val ans = f(a);map+=((a,ans));ans}
  }
}

class ParallelMemo[A,B](f:A=>B) extends Memo[A,B]{
  import scala.actors.Futures
  import scala.actors.Future
  var map=Map[A,Future[B]]()
  object Sync
  def future(a:A):Future[B]={
    if (map contains a){map(a)}
    else {
      Sync.synchronized{
        if (map contains a){map(a)} //check again in case of race condition
        else {
          val task = Futures.future{f(a)}
          map = map updated (a,task)
          task
        }
      }
    }
  }
  def apply(a:A):B={
    future(a).apply()
  }
  
}

class CacheMap[A,B](defaultSize:Int) extends scala.collection.mutable.Map[A,B]{
  val map = new org.apache.commons.collections.LRUMap(defaultSize)
    def get(a:A)={   
    val b = map.get(a)
      if (b==null){  
      None
    }else{
      Some(b.asInstanceOf[B])
    }
  }
  def -=(a:A)={
    map remove a
    this
  }
  def +=(kv:(A,B))={
    map put (kv._1,kv._2)
    this
  }
  def iterator={new JclIterator(map.iterator).map{_.asInstanceOf[A]}.map{a=>(a,get(a).get)}}
}
class JclIterator[A](i:java.util.Iterator[A]) extends Iterator[A]{
  def next=i.next
  def hasNext = i.hasNext
}


object Memo{
  def apply[A,B](f:A=>B)=new PermMemo(f)
}
object LRUMemo{
  def apply[A,B](size:Int)(f:A=>B)=new LRUMemo(size,f)
}
