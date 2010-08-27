package modiphy.util

trait Memo[A,B] extends Function[A,B]{
  def apply(a:A):B
}
class PermMemo[A, B](f:A=>B) extends Memo[A,B]{
  var map=Map[A,B]()
  def apply(a:A)={
    if (map contains a){map(a)}
    else {val ans = f(a);map=map updated (a,ans);ans}
  }
}
class LRUMemo[A,B](size:Int,f:A=>B) extends Memo[A,B]{
  val map=new CacheMap[A,B](size)
  def apply(a:A)={
    if (map contains a){map(a)}
    else {val ans = f(a);map+=((a,ans));ans}
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
