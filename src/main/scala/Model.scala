package modiphy.model
import modiphy.math.EnhancedMatrix._
import modiphy.tree._
import modiphy.util.Memo

trait Model{
  type Matrix=IndexedSeq[IndexedSeq[Double]]
  type Vector=IndexedSeq[Double]
  def apply(e:Edge):IndexedSeq[IndexedSeq[Double]]
  def pi(n:Node):IndexedSeq[Double]
}
trait StdModel extends Model{
  lazy val apply:Memo[Double,Matrix]=new Memo[Double,Matrix]({bl=>
    exp(bl)
  })
  def apply(e:Edge)=apply(e.dist)
  def exp(bl:Double):Matrix 
  def mat:Matrix
  /*
  def u:Matrix
  def d:Vector
  def v:Matrix
  */
}
trait ColtModel extends StdModel{
  import cern.colt.matrix.linalg._
  import cern.colt.matrix._
  import cern.colt.function._
  val algebra = new Algebra
  def expVals(d:DoubleMatrix2D,t:Double)=sparse.diagonal(fact2D.diagonal(d).assign( new DoubleFunction(){def apply(arg:Double)={math.exp(t * arg)}}))

  val eigen = try {new EigenvalueDecomposition(mat)}catch{ case e=> println("Can't find Eigensystem for Matrix " + mat + "\n because : " + e); throw e}
  val u = eigen.getV
  val v = algebra.inverse(u)
  val d = eigen.getD
  println("Mat " +mat)

  def exp(bl:Double):Matrix= { val ans = algebra.mult(algebra.mult(u,expVals(d,bl)),v);println("E " + bl + " " + ans);ans} 
}
abstract class ExpFactory{
  def apply(mat:Matrix)=make(mat)
  def make(mat:Matrix):Exp
}
abstract class CachingExpFactory extends ExpFactory{
  val cache:Function[Matrix,Exp]=LRUMemo[Matrix,Exp](50)({mat=>
    make(mat) 
  })
  override def apply(mat:Matrix)=cache(mat)
}
object ColtExpFactory extends CachingExpFactory{
  def make(mat:Matrix)=new ColtExp(mat)
}
object DefaultExpFactory extends ExpFactory{
  var default:ExpFactory=ColtExpFactory
  def make(mat:Matrix)=default.apply(mat)
  def setDefault(fact:ExpFactory){default=fact}
}

class BasicLikelihoodModel(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]],rate:Double=1.0,fact:ExpFactory=DefaultExpFactory) extends StdModel{
  val mat = s.sToQ(piValues,rate)
  def pi(n:Node) = piValues
  val exp=DefaultExpFactory.apply(mat)
}
>>>>>>> 3e29c6d... More debugging

}

abstract class BasicLikelihoodModel(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]]) extends StdModel{
  val mat = s.sToQ(piValues)
  def pi(n:Node) = piValues
}
