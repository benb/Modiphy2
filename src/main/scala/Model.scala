package modiphy.model
import modiphy.math.EnhancedMatrix._
import modiphy.tree._
import modiphy.util._
import scala.collection.LinearSeq
import modiphy.opt._


object Types{
  type Matrix=IndexedSeq[IndexedSeq[Double]]
  type LinearMatrix = LinearSeq[LinearSeq[Double]]
}
import Types._

trait Model{
  def paramMap:Map[(ParamName,Int),Any]
  lazy val params=paramMap.keys.toList
  def updatedVec(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]):Model
  def updatedMat(p:ParamName,vec:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int]):Model
  def getOptParam(p:ParamName,index:Option[Int]=None):IndexedSeq[Double]
  def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]):Model
}
trait SingleModel extends Model{
  val subModelIndex=0
  def pi(n:Node):IndexedSeq[Double]
  def apply(e:Edge):LinearMatrix
}
trait StdModel extends SingleModel{
  def pi:IndexedSeq[Double]
  def mat:Matrix
  val exp:Exp
  val rate=1.0
  def apply(e:Edge)=exp(e.dist * rate)
  lazy val paramMap=Map[(ParamName,Int),Any]((Pi,subModelIndex)->pi,(S,subModelIndex)->mat)
  def updatedVec(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]):StdModel
  def updatedMat(p:ParamName,vec:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int]):StdModel
  def updatedRate(d:Double):StdModel
  def getOptParam(p:ParamName,index:Option[Int]):IndexedSeq[Double]
  def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]):StdModel
  /*
  def u:Matrix
  def d:Vector
  def v:Matrix
  */
}
trait StdMixtureModel extends Model{
  def models:Seq[StdModel]
  def apply(e:Edge)=models.map{m=>m(e)}
  def priors:IndexedSeq[Double]
  lazy val paramMap = models.map{_.paramMap}.foldLeft(Map[(ParamName,Int),Any]()){_++_}
  def updatedVec(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]):StdMixtureModel
  def updatedMat(p:ParamName,vec:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int]):StdMixtureModel
  def getOptParam(p:ParamName,index:Option[Int]):IndexedSeq[Double]
  def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]):StdMixtureModel
}
trait Exp{
  def apply(bl:Double,rate:Double=1.0)=expInt(bl*rate)
  lazy val expInt:Memo[Double,LinearMatrix]=Memo[Double,LinearSeq[LinearSeq[Double]]]({bl=>
    exp(bl)
  })
  def exp(bl:Double):LinearSeq[LinearSeq[Double]]
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

  def exp(bl:Double):LinearMatrix=algebra.mult(algebra.mult(u,expVals(d,bl)),v) 
}
abstract class ExpFactory{
  def apply(mat:Matrix):Exp=make(mat)
  def apply(mat:Matrix,pi:IndexedSeq[Double]):Exp=apply(mat.sToQ(pi))
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

object BasicLikelihoodModel{
  def apply(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]]):BasicLikelihoodModel={
    apply(piValues,s,0)
  }
  def apply(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]],subModelIndex:Int):BasicLikelihoodModel={
    apply(piValues,s,1.0,subModelIndex)
  }
  def apply(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]],rate:Double,subModelIndex:Int):BasicLikelihoodModel={
    new BasicLikelihoodModel(piValues,s,rate,DefaultExpFactory(s,piValues),subModelIndex)
  }
}

class BasicLikelihoodModel(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]], override val rate:Double, val exp:Exp,override val subModelIndex:Int=0) extends StdModel{
  val mat = s.sToQ(piValues,rate)
  def pi(n:Node) = piValues
  val pi = piValues
  def updatedRate(r:Double)={
     new BasicLikelihoodModel(piValues,s,r,exp,subModelIndex)
  }
  def updatedVec(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int])={
    val mI = paramIndex getOrElse subModelIndex
    (p,mI) match {
      case (Pi,`subModelIndex`) => BasicLikelihoodModel(vec,s,rate,subModelIndex) 
    }
  }
  def updatedMat(p:ParamName,vec:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int])= {
    val mI = paramIndex getOrElse subModelIndex
    (p,mI) match {
      case (S,`subModelIndex`) => BasicLikelihoodModel(piValues,vec,rate,subModelIndex)
    }
  }
  def getOptParam(p:ParamName,modelIndex:Option[Int]):IndexedSeq[Double]={
    val mI = modelIndex.getOrElse(subModelIndex)
    (p,mI) match {
      case (Pi,`subModelIndex`)=> Pi.getOpt(piValues)
      case (S,`subModelIndex`) => S.getOpt(s)
    }
  }
  def setOptParam(p:ParamName,vec:IndexedSeq[Double],modelIndex:Option[Int])={
    val mI = modelIndex.getOrElse(subModelIndex)
    (p,mI) match {
      case (Pi,`subModelIndex`)=> BasicLikelihoodModel(Pi.getReal(vec),s,rate,subModelIndex)
      case (S,`subModelIndex`)=> BasicLikelihoodModel(piValues,S.getReal(vec),rate,subModelIndex)
    }
  }
}

object GammaModel{
  def apply(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]],alpha:Double,numCat:Int):GammaModel={
    apply(piValues,s,alpha,numCat,0)
  }
  def apply(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]],alpha:Double,numCat:Int,myParamIndex:Int):GammaModel={
    apply(piValues,s,alpha,numCat,myParamIndex,BasicLikelihoodModel(piValues,s))
  }
  def apply(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]],alpha:Double,numCat:Int,myParamIndex:Int,basicModels1:StdModel):GammaModel={
    apply(piValues,s,alpha,numCat,myParamIndex,Gamma.getDist(alpha,numCat).map{t=>basicModels1.updatedRate(t)})  
  }
  def apply(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]],alpha:Double,numCat:Int,myParamIndex:Int,models:Seq[StdModel]):GammaModel={
    new GammaModel(piValues,s,alpha,numCat,myParamIndex,models)
  }
}
class GammaModel(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]],alpha:Double,numCat:Int,myParamIndex:Int,val models:Seq[StdModel]) extends StdMixtureModel{
  val priors = Vector.fill(numCat)(1.0/numCat)
  //stub methods FIXME
  def updatedVec(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]) = {
    val mI = paramIndex.getOrElse(myParamIndex)
    (p,mI) match {
      case (Gamma,`myParamIndex`)=>GammaModel(piValues,s,vec(0),numCat,myParamIndex,models.head)
      case (Pi,`myParamIndex`)=>GammaModel(vec,s,alpha,numCat,myParamIndex)
      case _ => this
    }
  }
  def updatedMat(p:ParamName,vec:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int]) = {
    val mI = paramIndex.getOrElse(myParamIndex)
    (p,mI) match {
      case (S,`myParamIndex`)=>GammaModel(piValues,vec,alpha,numCat,myParamIndex)
      case _ => this
    }
  }

  def getOptParam(p:ParamName,paramIndex:Option[Int]):IndexedSeq[Double]={
    val mI = paramIndex.getOrElse(myParamIndex)
    (p,mI) match {
      case (Pi,`myParamIndex`)=> Pi.getOpt(piValues)
      case (S,`myParamIndex`)=> S.getOpt(s)
      case (Gamma,`myParamIndex`) => Vector(alpha)
    }
  }
  def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int])={
    val mI = paramIndex.getOrElse(myParamIndex)
    (p,mI) match {
      case (Pi,`myParamIndex`)=> updatedVec(Pi,Pi.getReal(vec),paramIndex)
      case (S,`myParamIndex`)=> updatedMat(S,S.getReal(vec),paramIndex)
      case (Gamma,`myParamIndex`) => updatedVec(Gamma,vec,paramIndex)
      case _ => this
    }
 
  }

}

class GammaDist{
  import org.apache.commons.math.distribution.ChiSquaredDistributionImpl
  import cern.jet.stat.Gamma.incompleteGamma
  val chi2=new ChiSquaredDistributionImpl(1.0D)
  def chiSquareInverseCDF(prob:Double,df:Double)={
     chi2.setDegreesOfFreedom(df)
     val ans = chi2.inverseCumulativeProbability(prob)
     ans
  }
  def gammaInverseCDF(prob:Double,alpha:Double,beta:Double)=chiSquareInverseCDF(prob,2.0*(alpha))/(2.0*(beta))


  def apply(shape:Double,numCat:Int):IndexedSeq[Double]={
    gamma(shape,numCat)
  }
  def gamma(shape:Double,numCat:Int):IndexedSeq[Double]={
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

