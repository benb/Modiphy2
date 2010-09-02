package modiphy.model
import modiphy.math.EnhancedMatrix._
import modiphy.tree._
import modiphy.util._
import scala.collection.LinearSeq
import modiphy.opt._


object Types{
  type Matrix=LinearSeq[LinearSeq[Double]]
}
import Types._

trait Model{
  def paramMap:Map[(ParamName,Int),Any]
  lazy val params=paramMap.keys.toList
}
trait SingleModel extends Model{
  val subModelIndex=0
  def pi(n:Node):IndexedSeq[Double]
  def apply(e:Edge):Matrix
}
trait StdModel extends SingleModel{
  def pi:IndexedSeq[Double]
  def mat:Matrix
  val exp:Exp
  def apply(e:Edge)=exp(e.dist)
  lazy val paramMap=Map[(ParamName,Int),Any]((Pi,subModelIndex)->pi,(S,subModelIndex)->mat)
  /*
  def u:Matrix
  def d:Vector
  def v:Matrix
  */
}
trait StdMixtureModel extends Model{
  def models:Seq[SingleModel]
  def apply(e:Edge)=models.map{m=>m(e)}
  def priors:IndexedSeq[Double]
  lazy val paramMap = models.map{_.paramMap}.foldLeft(Map[(ParamName,Int),Any]()){_++_}
}
trait Exp{
  def apply(bl:Double,rate:Double=1.0)=expInt(bl*rate)
  lazy val expInt:Memo[Double,Matrix]=Memo[Double,Matrix]({bl=>
    exp(bl)
  })
  def exp(bl:Double):Matrix
}

class ColtExp(mat:Matrix) extends Exp{
  import cern.colt.matrix.linalg._
  import cern.colt.matrix._
  import cern.colt.function._
  val algebra = new Algebra
  def expVals(d:DoubleMatrix2D,t:Double)=sparse.diagonal(fact2D.diagonal(d).assign( new DoubleFunction(){def apply(arg:Double)={math.exp(t * arg)}}))

  val eigen = try {new EigenvalueDecomposition(mat)}catch{ case e=> println("Can't find Eigensystem for Matrix " + mat + "\n because : " + e); throw e}
  val u = eigen.getV
  val v = algebra.inverse(u)
  val d = eigen.getD

  def exp(bl:Double):Matrix=algebra.mult(algebra.mult(u,expVals(d,bl)),v) 
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

class BasicLikelihoodModel(piValues:IndexedSeq[Double],s:LinearSeq[LinearSeq[Double]],rate:Double=1.0,fact:ExpFactory=DefaultExpFactory,override val subModelIndex:Int=0) extends StdModel{
  val mat = s.sToQ(piValues,rate)
  def pi(n:Node) = piValues
  val pi = piValues
  val exp=DefaultExpFactory.apply(mat)
  
}

class GammaModel(piValues:IndexedSeq[Double],s:LinearSeq[LinearSeq[Double]],alpha:Double,numCat:Int) extends StdMixtureModel{
  lazy val gamma = new Gamma(numCat).apply(alpha)
  val models:Seq[StdModel] = gamma.zipWithIndex.map{t=>new BasicLikelihoodModel(piValues,s,t._1,subModelIndex=t._2)}
  val priors = Vector.fill(numCat)(1.0/numCat)
}

class Gamma(numCat:Int){
  import org.apache.commons.math.distribution.ChiSquaredDistributionImpl
  import cern.jet.stat.Gamma.incompleteGamma
  val chi2=new ChiSquaredDistributionImpl(1.0D)
  def chiSquareInverseCDF(prob:Double,df:Double)={
     chi2.setDegreesOfFreedom(df)
     val ans = chi2.inverseCumulativeProbability(prob)
     ans
  }
  def gammaInverseCDF(prob:Double,alpha:Double,beta:Double)=chiSquareInverseCDF(prob,2.0*(alpha))/(2.0*(beta))


  def apply(shape:Double):IndexedSeq[Double]={
    gamma(shape)
  }
  def gamma(shape:Double):IndexedSeq[Double]={
    if (shape==Double.PositiveInfinity){
      (0 until numCat).map{a=>1.0}
    }else {
      val alpha = shape
      val beta = shape
      val factor=alpha/beta*numCat
      val freqK=new Array[Double](numCat)
        val rK=new Array[Double](numCat)

        (0 until numCat-1).foreach{i=>
          freqK(i)=gammaInverseCDF((i+1.0)/numCat, alpha, beta);
          freqK(i) = incompleteGamma(alpha+1,freqK(i)*beta)
        }

        rK(0) = freqK(0)*factor;
        rK(numCat-1) = (1-freqK(numCat-2))*factor;
        (1 until numCat-1).foreach{i=>
        rK(i) = (freqK(i)-freqK(i-1))*factor;
      }
      // println("RATES " + rK.toList)
      rK.toIndexedSeq
    }
  }
}

