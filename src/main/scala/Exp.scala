package modiphy.model
import modiphy.math.EnhancedMatrix._
import modiphy.tree._
import modiphy.util._
import scala.collection.LinearSeq
import modiphy.opt._
import modiphy.alignment.Alphabet
import modiphy.math._

import Types._
trait Exp{
  def realExp = this
  def apply(bl:Double,rate:Double=1.0)=expInt(bl*rate)
  lazy val expInt:Memo[Double,LinearSeq[LinearSeq[Double]]]=Memo[Double,LinearSeq[LinearSeq[Double]]]({bl=>
    exp(bl)
  })
  def exp(bl:Double):LinearSeq[LinearSeq[Double]]
  def rerate(newRate:Double):Exp = {
    val outer = realExp
    new Exp{
      override def realExp = outer
      def exp(bl:Double)=outer.exp(bl*newRate)
    }
  }
  def instantiated:Boolean=true
}

class ColtExp(mat:Matrix,rate:Double=1.0) extends Exp{
  import cern.colt.matrix.linalg._
  import cern.colt.matrix._
  import cern.colt.function._
  import Timing.time
  lazy val algebra = new Algebra
  def expVals(d:DoubleMatrix2D,t:Double)=sparse.diagonal(fact2D.diagonal(d).assign( new DoubleFunction(){def apply(arg:Double)={math.exp(t * arg)}}))

  lazy val eigen = {
    val ans = time{
      try {new EigenvalueDecomposition(mat)}catch{ case e=> println("Can't find Eigensystem for Matrix " + mat + "\n because : " + e); throw e}
    }
    ans._1
  }
  lazy val u = eigen.getV
  lazy val v = algebra.inverse(u)
  lazy val d = eigen.getD

  def exp(bl:Double):LinearMatrix=algebra.mult(algebra.mult(u,expVals(d,bl)),v) 
  override def toString = {
    "Left Eigenvectors " + u + "\nEigenvalues" + d  + "\nRight Eigenvectors " + v
  }
}
class JBlasExp(m:Matrix) extends Exp{
  import org.jblas.{DoubleMatrix,Solve,MatrixFunctions}
  import org.jblas.Eigen._
  import Timing.time

  lazy val mat = new DoubleMatrix(m.length,m.length,m.flatten:_*)
  lazy val eigen = {
    val ans = time{
      eigenvectors(mat)
    }
    ans._1
  }
  lazy val u = eigen(0).getReal
  lazy val d = eigen(1).getReal
  lazy val v = Solve.solve(u.dup,DoubleMatrix.eye(u.getRows))
  lazy val size = mat.getRows

  def exp(bl:Double):LinearMatrix = {
    val myD = DoubleMatrix.diag(MatrixFunctions.expi(d.mmul(bl)).diag)
    val ans = u.mmul(myD).mmul(v).transpose
    val out = ans.toArray2.map{_.toList}.toList
    out
  }

}

object Timing{
  def time[A](f: =>A):(A,Double)={
    val start = System.nanoTime : Double
    val ans = f
    val end = System.nanoTime : Double
    (ans,(end-start)/1000000.0)
  }
}

abstract class ExpFactory{
  def apply(mat:Matrix):Exp=make(mat)
  def apply(mat:Matrix,pi:IndexedSeq[Double],rate:Double):Exp=apply(mat.sToQ(pi,rate))
  def make(mat:Matrix):Exp
}
abstract class CachingExpFactory extends ExpFactory{
  val cache:Function[Matrix,Exp]=LRUMemo[Matrix,Exp](50)({mat=>
    make(mat) 
  })
  override def apply(mat:Matrix)=cache(mat)
}
object ColtExpFactory extends CachingExpFactory{
  def make(mat:Matrix):ColtExp=new ColtExp(mat)
}
object JBlasExpFactory extends CachingExpFactory{
  def make(mat:Matrix)=new JBlasExp(mat)
}
object DefaultExpFactory extends ExpFactory{
  var default:ExpFactory=ColtExpFactory
  def make(mat:Matrix)=default.apply(mat)
  def setDefault(fact:ExpFactory){default=fact}
}


