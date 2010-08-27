package modiphy.math
import scala.collection.LinearSeq
class EnhancedMatrix(mat:IndexedSeq[IndexedSeq[Double]]){
  implicit def MakeEnhancedMatrix(mat:IndexedSeq[IndexedSeq[Double]])=new EnhancedMatrix(mat)
  def sToQ(pi:IndexedSeq[Double],rate:Double=1.0)={
    IndexedSeq.tabulate(mat.length,mat.head.length){(i,j)=>
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
}
class EnhancedVector(seq:IndexedSeq[Double]){
  def normalize = normalise
  def normalise = {
    val mySum = sum
    seq.map{i=> i/sum}
  }
  def sum = seq.reduceLeft{_+_}
  def dotProduct(vect:IndexedSeq[Double])=seq.zip(vect).map{t=>t._1*t._2}.reduceLeft{_+_}
}
class EnhancedListVector(seq:LinearSeq[Double]){
  def normalize = normalise
  def normalise = {
    val mySum = sum
    seq.map{i=> i/sum}
  }
  def sum = seq.reduceLeft{_+_}
  def dotProduct(vect:IndexedSeq[Double])=seq.zip(vect).map{t=>t._1*t._2}.reduceLeft{_+_}
}

object EnhancedMatrix{
  import cern.colt.matrix._
  implicit def MakeEnhancedMatrix(mat:IndexedSeq[IndexedSeq[Double]])=new EnhancedMatrix(mat)
  implicit def MakeEnhancedVector(vec:IndexedSeq[Double])=new EnhancedVector(vec)
  implicit def MakeEnhancedListVector(vec:LinearSeq[Double])=new EnhancedListVector(vec)
  lazy val fact1D = cern.colt.matrix.DoubleFactory1D.dense
  lazy val fact2D = cern.colt.matrix.DoubleFactory2D.dense
  lazy val sparse = cern.colt.matrix.DoubleFactory2D.sparse

  //lazy initial pass at these methods
  implicit def Seq2Vec(s:IndexedSeq[Double])=fact1D.make(s.toArray)
  implicit def Vec2Seq(v:DoubleMatrix1D)=v.toArray.toIndexedSeq
  implicit def Seq2Mat(s:IndexedSeq[IndexedSeq[Double]])=fact2D.make(s.map{_.toArray}.toArray)
  implicit def Mat2Seq(m:DoubleMatrix2D)=m.toArray.map{_.toIndexedSeq}.toIndexedSeq
}
