package modiphy.alignment
import scala.collection.immutable.VectorBuilder

abstract class Alphabet(names:String*) extends Enumeration(names :_*) {
  import GlobalAlphabet._
  class Letter(name:String,val alphabet:Alphabet) extends Val(nextId,name)
  def Letter = new Letter(if (nextName.hasNext) nextName.next else null,this)
  def unknown:Letter
  def isReal(v:Letter):Boolean
  def matElements:Seq[Letter]
  val length = matElements.length
  def gap:Letter
  def parseString(s:String):IndexedSeq[Letter]
}


object GlobalAlphabet{
  type Letter=Alphabet#Letter
  case class EnhancedLetter[A<:Alphabet](l:A#Letter){
    def isReal=l.alphabet.isReal(l.asInstanceOf[l.alphabet.Letter])
  }
  implicit def MakeLetter(l:Alphabet#Letter)=EnhancedLetter(l)
}

import GlobalAlphabet._


object DNA extends Alphabet("A","G","C","T","N","-"){
  val A,G,C,T,N,GAP=Letter
  def isReal(v:Letter)={v!=N && v!=GAP}
  def parseString(s:String)=s.map{i=> values.find(v=>v.toString==i.toString).getOrElse(N)}.foldLeft(new VectorBuilder[Letter]){_ += _.asInstanceOf[Letter]}.result
  def matElements=List(A,G,C,T)
  val unknown=N
  val gap=GAP
}

object AminoAcid extends Alphabet("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","X","-"){
  val A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,X,GAP=Letter
  override def isReal(a:Letter)=((a!=X) && (a!=GAP))
  def matElements=List(A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V)
  def parseString(s:String)=s.map{i=> values.find(v=>v.toString==i.toString).getOrElse(X)}.map{_.asInstanceOf[Letter]}.toIndexedSeq
  val unknown=X
  val gap=GAP
}


class Aligmment(gen:Map[String,Iterable[Letter]]){
  val names = gen.keys.toList
  type Pattern=Map[String,Letter]
  val columns:Seq[Pattern] = FlippedIterator(names.map{gen}.map{_.iterator}).map{_.toList}.foldLeft(List[Map[String,Letter]]()){(ml,col)=>
    names.zip(col).foldLeft(Map[String,Letter]()){(ml2,item)=>
      ml2 + item
    } :: ml
  }.reverse
  lazy val patterns:Map[Pattern,Int]=columns.foldLeft(Map[Pattern,Int]()){(m,col)=>
    m.updated(col,m.getOrElse(col,0)+1)
  }
  lazy val patternList=patterns.toList.map{_._1}
  lazy val countList=patterns.toList.map{_._2}
  def apply(s:String)=gen(s)
  lazy val frequencies={
    val alphabet = columns.head.values.head.alphabet
    var countMap=Map[Letter,Int]()
    columns.map{p:Pattern=>
      countMap = p.values.foldLeft(countMap){(m,letter)=>
        m updated (letter,m.getOrElse(letter,0)+1) 
      }
    }
    val total = countMap.filter{_._1.isReal}.values.reduceLeft(_+_)
    alphabet.matElements.map{countMap.getOrElse(_,0)}.map{_.toDouble/total}
  }.toIndexedSeq
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
  def parseWith(alphabet:Alphabet)={
    val map = this.foldLeft(Map[String,Seq[Letter]]()){(m,t)=>
      m updated (t._1,alphabet.parseString(t._2)) 
    }
    new Aligmment(map)
  }
}
object Fasta{
  def apply(source:Iterator[String])=new Fasta(source)
  def apply(source:String)=new Fasta(source.lines)
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

