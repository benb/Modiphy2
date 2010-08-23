package modiphy.alignment
import scala.collection.immutable.VectorBuilder

object GlobalAlphabet{
  case class EnhancedLetter[A<:Alphabet](l:A#Value,a:A){
    def id:Int=l.id
    def alphabet:Alphabet=a
    def isReal=a.isReal(l.asInstanceOf[a.Value])
  }

  private var letterMap=Map[Alphabet#Value,EnhancedLetter[_]]()
  def addNew(a:Alphabet){
    a.values.foreach{v=>
      letterMap = letterMap updated (v,EnhancedLetter(v,a))
    }
  }
  implicit def MakeLetter(l:Alphabet#Value)=letterMap(l)
  type Letter=Alphabet#Value
}

import GlobalAlphabet._

abstract class Alphabet(names:String*) extends Enumeration(names :_*) {
  import GlobalAlphabet._
  addNew(this)
  def unknown:Value
  def isReal(v:Value):Boolean
  def matElements:Seq[Value]
  val length = matElements.length
  def gap:Value
  def parseString(s:String):IndexedSeq[Value]
}

object DNA extends Alphabet("A","G","C","T","N","-"){
  val A,G,C,T,N,GAP=Value
  def isReal(v:Value)={v!=N && v!=GAP}
  def parseString(s:String)=s.map{i=> values.find(v=>v.toString==i.toString).getOrElse(N)}.foldLeft(new VectorBuilder[Value]){_ += _}.result
  def matElements=List(A,G,C,T)
  val unknown=N
  val gap=GAP
}

object AminoAcid extends Alphabet("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","X","-"){
  val A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,X,GAP=Value
  override def isReal(a:Value)=((a!=X) && (a!=GAP))
  def matElements=List(A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V)
  def parseString(s:String)=s.map{i=> values.find(v=>v.toString==i.toString).getOrElse(X)}.toIndexedSeq
  val unknown=X
  val gap=GAP
}


class Aligmment(gen:Map[String,Iterable[Letter]]){
  val names = gen.keys.toList
  val columns:Seq[Map[String,Letter]] = FlippedIterator(names.map{gen}.map{_.iterator}).map{_.toList}.foldLeft(List[Map[String,Letter]]()){(ml,col)=>
    names.zip(col).foldLeft(Map[String,Letter]()){(ml2,item)=>
      ml2 + item
    } :: ml
  }.reverse
  def apply(s:String)=gen(s)
}

class Fasta(source:Iterator[String]) extends Iterator[(String,String)]{
  import EnhancedIterator._
  val iter=source.map{i=>i.trim}.buffered
  def next = {
    val name = iter.next.trim.drop(1).trim.split("\\s+")(0)
    val seq = iter.takeWhileSafe{s:String => !(s startsWith (">"))}.mkString("").toUpperCase
    (name,seq)
  }
  def hasNext = iter.hasNext
  def toAlignment(alphabet:Alphabet)={
    val map = this.foldLeft(Map[String,Seq[Letter]]()){(m,t)=>
      m updated (t._1,alphabet.parseString(t._2)) 
    }
    new Aligmment(map)
  }
}


object FlippedIterator{
  def apply[A](l:Iterable[Iterator[A]])(implicit m:Manifest[Iterable[Iterator[A]]])=new FlippedIterator(l)
    implicit def MakeIterator[A](l:Iterable[Iterable[A]])=l.map{_.iterator}
  //  def apply[A](l:Iterable[Iterable[A]])(implicit m:Manifest[Iterable[Iterable[A]]])=new FlippedIterator(l.map{_.iterator})
}
class FlippedIterator[A](l:Iterable[Iterator[A]]) extends Iterator[Iterable[A]]{
  def next = l.map{_.next}
  def hasNext = l.foldLeft(true){(a,b)=>a && b.hasNext}
}

object EnhancedIterator{
  implicit def MakeEnhancedIterator[A](i: Iterator[A]):EnhancedIterator[A]=new EnhancedIterator(i.buffered)
}
class EnhancedIterator[A](i:BufferedIterator[A]){
  def takeWhileSafe(f:A=>Boolean):List[A]={
    takeWhileSafe(List(),f)
  }
  def takeWhileSafe(l:List[A],f:A=>Boolean):List[A] ={
    if (i.hasNext && f(i.head)){
      takeWhileSafe(i.next::l,f)
    }else {
      l.reverse
    }
  }
}

