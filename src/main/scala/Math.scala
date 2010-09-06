package modiphy.math
import scala.collection.LinearSeq
import EnhancedMatrix._
class EnhancedIndexedMatrix(mat:IndexedSeq[IndexedSeq[Double]]){
  def sToQ(pi:IndexedSeq[Double],rate:Double=1.0)={
    //fixme
   List.tabulate(mat.length,mat.head.length){(i,j)=>
      if (j>i){
        mat(j)(i)*pi(j)
      }else{
        mat(i)(j)*pi(j)
      }
    }.fixDiag.normalise(pi,rate)
  }
  def fixDiag={
    mat.zipWithIndex.map{t=> val (row,i)=t
      val sum= row.updated(i,0.0).reduceLeft{_+_}
      row.updated(i,-sum)
    }
  }
  def normalise(pi:IndexedSeq[Double],rate:Double=1.0)={
    val zSum = mat.diag.zip(pi).map{t=>t._1*t._2}.sum
    mat.map{row=>
      row.map{_ / -zSum * rate}
    }
  }
  def diag={
    mat.zipWithIndex.map{t=> val (row,i)=t
      row(i)
    }
  }
}
class EnhancedMatrix(mat:LinearSeq[LinearSeq[Double]]){
  import EnhancedMatrix._
  def sToQ(pi:IndexedSeq[Double],rate:Double=1.0)={
    //fixme
    List.tabulate(mat.length,mat.head.length){(i,j)=>
      if (j>i){
        mat(j)(i)*pi(j)
      }else{
        mat(i)(j)*pi(j)
      }
    }.fixDiag.normalise(pi,rate)
  }
  def fixDiag={
    mat.zipWithIndex.map{t=> val (row,i)=t
      val sum= row.updated(i,0.0).reduceLeft{_+_}
      row.updated(i,-sum)
    }
  }
  def diag={
    mat.zipWithIndex.map{t=> val (row,i)=t
      row(i)
    }
  }
  def normalise(pi:IndexedSeq[Double],rate:Double=1.0)={
    val zSum = mat.diag.zip(pi).map{t=>t._1*t._2}.sum
    mat.map{row=>
      row.map{_ / -zSum * rate}
    }
  }
  def multipleDotProduct(mat2:LinearSeq[LinearSeq[Double]]):LinearSeq[LinearSeq[Double]] = {
    mat.map{e=>
      mat2.map{e.dotProduct}
    }
  }

}
trait EnhancedVector[A <: Seq[Double]]{
  def normalise:A
  def normalize:A = normalise
  def sum:Double
  def dotProduct(other:Seq[Double]):Double
  def prod(other:Seq[Double]):A
  
}
class EnhancedIndexedVector(seq:IndexedSeq[Double]) extends EnhancedVector[IndexedSeq[Double]]{
  def normalise = {
    val mySum = sum
    seq.map{i=> i/sum}
  }
  def sum = seq.reduceLeft{_+_}
  def dotProduct(vect:Seq[Double])={
    var ans = 0.0
    seq.zip(vect).foreach{t=>ans = ans + t._1*t._2}
    ans
  }
  def prod(vect:Seq[Double])={
    seq.zip(vect).map{t=>t._1*t._2}
  }
}
class EnhancedListVector(seq:LinearSeq[Double]) extends EnhancedVector[LinearSeq[Double]]{
  def normalise = {
    val mySum = sum
    seq.map{i=> i/sum}
  }
  def sum = seq.reduceLeft{_+_}
  def dotProduct(vect:Seq[Double])={
    var ans = vect.head * seq.head
    var l1 = vect.tail
    var l2 = seq.tail
    while (!(l1.isEmpty)){
      ans = ans + l1.head * l2.head
      l1 = l1.tail
      l2 = l2.tail
    }
    ans
  }
  def prod(vect:Seq[Double])={
    var ans = vect.head * seq.head :: Nil
    var l1 = vect.tail
    var l2 = seq.tail
    while (!(l1.isEmpty)){
      ans = l1.head * l2.head :: ans
      l1 = l1.tail
      l2 = l2.tail
    }
    ans.reverse
  }
}

object EnhancedMatrix{
  import cern.colt.matrix._
  implicit def MakeEnhancedMatrix(mat:LinearSeq[LinearSeq[Double]])=new EnhancedMatrix(mat)
  implicit def MakeEnhancedIndexedMatrix(mat:IndexedSeq[IndexedSeq[Double]])=new EnhancedIndexedMatrix(mat)
  implicit def MakeEnhancedListVector(vec:LinearSeq[Double]) = new EnhancedListVector(vec)
  implicit def MakeEnhancedIndexedVector(vec:IndexedSeq[Double]) = new EnhancedIndexedVector(vec)

  lazy val fact1D = cern.colt.matrix.DoubleFactory1D.dense
  lazy val fact2D = cern.colt.matrix.DoubleFactory2D.dense
  lazy val sparse = cern.colt.matrix.DoubleFactory2D.sparse

  //lazy initial pass at these methods
  implicit def Seq2Vec(s:Seq[Double])=fact1D.make(s.toArray)
  implicit def Vec2Seq(v:DoubleMatrix1D)=v.toArray.toList
 // implicit def Vec2LinearSeq(v:DoubleMatrix1D)=v.toArray.toList
  implicit def Seq2Mat(s:LinearSeq[LinearSeq[Double]])=fact2D.make(s.map{_.toArray}.toArray)
  implicit def Mat2Seq(m:DoubleMatrix2D)=m.toArray.map{_.toList}.toList
}
