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
  def getOptParam(p:ParamName,index:Option[Int]=None):Option[IndexedSeq[Double]]
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
  def getOptParam(p:ParamName,index:Option[Int]):Option[IndexedSeq[Double]]
  def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]):StdModel
  def numClasses = 1
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
  lazy val paramMap = models.map{_.paramMap}.foldLeft(mixtureParams){_++_}
  def mixtureParams=Map[(ParamName,Int),Any]()
  def updatedVec(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]):StdMixtureModel
  def updatedMat(p:ParamName,vec:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int]):StdMixtureModel
  def getOptParam(p:ParamName,index:Option[Int]):Option[IndexedSeq[Double]]
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
  lazy val algebra = new Algebra
  def expVals(d:DoubleMatrix2D,t:Double)=sparse.diagonal(fact2D.diagonal(d).assign( new DoubleFunction(){def apply(arg:Double)={math.exp(t * arg)}}))

  lazy val eigen = try {new EigenvalueDecomposition(mat)}catch{ case e=> println("Can't find Eigensystem for Matrix " + mat + "\n because : " + e); throw e}
  lazy val u = eigen.getV
  lazy val v = algebra.inverse(u)
  lazy val d = eigen.getD

  def exp(bl:Double):LinearMatrix=algebra.mult(algebra.mult(u,expVals(d,bl)),v) 
}
class JBlasExp(m:Matrix) extends Exp{
  import org.jblas.{DoubleMatrix,Solve,MatrixFunctions}
  import org.jblas.Eigen._

  lazy val mat = new DoubleMatrix(m.length,m.length,m.flatten:_*)
  lazy val eigen = symmetricEigenvectors(mat)
  lazy val u = eigen(0)
  lazy val d = eigen(1)
  lazy val v = Solve.solveSymmetric(u,DoubleMatrix.eye(u.getRows))

  def exp(bl:Double):LinearMatrix = {
    val ans = u.dup
    val myD = MatrixFunctions.expi(d.dup.mmul(bl))
    u.mmul(myD).mmul(v)
    ans.toArray2.toList.map{_.toList}
  }

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
object JBlasExpFactory extends CachingExpFactory{
  def make(mat:Matrix)=new JBlasExp(mat)
}
object DefaultExpFactory extends ExpFactory{
  var default:ExpFactory=JBlasExpFactory
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
/*
  Site class model with one n*m n*m matrix for n site classes of m-length alphabet.
  Subclasses could implement Covarion or THMM type models.
*/
object StdSiteClassModel{
  def apply(subModel:Seq[StdModel],priors:IndexedSeq[Double],modelIndex:Int,bigMat:IndexedSeq[IndexedSeq[Double]]) = {
    new StdSiteClassModel(subModel,priors,modelIndex,Some(bigMat),None)
  }
  def apply(subModel:Seq[StdModel],priors:IndexedSeq[Double],modelIndex:Int)={
    new StdSiteClassModel(subModel,priors,modelIndex,None,None)
  }
  def apply(gamma:GammaModel):StdSiteClassModel={ apply(gamma.models) }
  def apply(models:Seq[StdModel]):StdSiteClassModel={
    new StdSiteClassModel(models,Vector.fill(models.length)(1.0/models.length),0,None,None)
  }
}
class StdSiteClassModel(subModel:Seq[StdModel],priors:IndexedSeq[Double],modelIndex:Int=0,var bigMat:Option[IndexedSeq[IndexedSeq[Double]]],var myExp:Option[Exp]) extends StdModel{
  override val numClasses = priors.length
  lazy val mat = bigMat getOrElse {
    val ans = subModel.map{_.mat}.foldLeft[IndexedSeq[IndexedSeq[Double]]]( Vector[Vector[Double]]() ){(m,m2)=>m addClass m2}
    bigMat = Some(ans)
    ans
  }
  lazy val pi = subModel.map{_.pi}.zip(priors).map{t=> t._1.map{_*t._2}}.flatten.toIndexedSeq
  lazy val exp = myExp.getOrElse{val ans = DefaultExpFactory(mat); myExp = Some(ans); ans}
  def pi(n:Node)=pi
  def updatedVec(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int])={
    p match {
      case MixturePrior if (paramIndex.isEmpty || paramIndex.get==modelIndex) => new StdSiteClassModel(subModel,vec,modelIndex,bigMat,myExp)
      case any => new StdSiteClassModel(subModel.map{_.updatedVec(p,vec,paramIndex)},priors,modelIndex,None,None)
    }
  }
  def updatedMat(p:ParamName,vec:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int])= {
     new StdSiteClassModel(subModel.map{_.updatedMat(p,vec,paramIndex)},priors,modelIndex,bigMat,myExp)
   }
   def getOptParam(p:ParamName,paramIndex:Option[Int])={
     if (p==MixturePrior && (paramIndex.isEmpty || paramIndex.get==modelIndex)){
       Some(MixturePrior.getOpt(priors))
     }else {
       subModel.foldLeft[Option[IndexedSeq[Double]]](None){(ans,model)=>
         if (ans.isEmpty){
           model.getOptParam(p,paramIndex)
         }else {
           ans
         }
       }
     }
   }
  def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int])={
    p match {
      case MixturePrior if (paramIndex.isEmpty || paramIndex.get==modelIndex) => new StdSiteClassModel(subModel,MixturePrior.getReal(vec),modelIndex,bigMat,myExp)
      case any => new StdSiteClassModel(subModel.map{_.setOptParam(p,vec,paramIndex)},priors,modelIndex,None,None)
    }
  }
  def updatedRate(r:Double):StdSiteClassModel={error ("Rate currently unimplemented for StdSiteClassModel"); this}
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
  def getOptParam(p:ParamName,modelIndex:Option[Int])={
    val mI = modelIndex.getOrElse(subModelIndex)
    (p,mI) match {
      case (Pi,`subModelIndex`)=> Some(Pi.getOpt(piValues))
      case (S,`subModelIndex`) => Some(S.getOpt(s))
      case (_,_) => None
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
  override def mixtureParams=Map((Gamma,myParamIndex)->alpha)
  //stub methods FIXME
  def updatedVec(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]) = {
    val mI = paramIndex.getOrElse(myParamIndex)
    (p,mI) match {
      case (Gamma,`myParamIndex`)=>GammaModel(piValues,s,vec(0),numCat,myParamIndex,models.head)
      case (Pi,`myParamIndex`)=>{
        GammaModel(vec,s,alpha,numCat,myParamIndex)
      }
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

  def getOptParam(p:ParamName,paramIndex:Option[Int])={
    val mI = paramIndex.getOrElse(myParamIndex)
    (p,mI) match {
      case (Pi,`myParamIndex`)=> Some(Pi.getOpt(piValues))
      case (S,`myParamIndex`)=> Some(S.getOpt(s))
      case (Gamma,`myParamIndex`) => Some(Vector(alpha))
      case (_,_) => None
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

