package modiphy.model
import modiphy.math.EnhancedMatrix._
import modiphy.tree._
import modiphy.util._
import scala.collection.LinearSeq
import modiphy.opt._
import modiphy.alignment.Alphabet
import modiphy.math._

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
  def updatedRate(r:Double):Model
  def numClasses:Int
  def mats:Seq[Matrix]
  def pis:Seq[IndexedSeq[Double]]
  def pi:IndexedSeq[Double]
  def mat:Matrix
  def pi(n:Node):IndexedSeq[Double]=pi
  val exp:Exp
  def priors:IndexedSeq[Double]
  val rate=1.0
  def apply(e:Edge)=exp(e.dist * rate)
  def models:Seq[Model]=List(this)

}
trait SingleModel extends Model{
  val subModelIndex=0
  def pi(n:Node):IndexedSeq[Double]
  def apply(e:Edge):LinearMatrix
}
trait StdModel extends SingleModel{
  def priors=Vector(1.0)
  def pi:IndexedSeq[Double]
  def mat:Matrix
  def mats=List(mat)
  def pis = List(pi)
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
  override def models:Seq[Model]
  def mats=models.map{_.mats}.flatten
  def pis = models.map{_.pis}.flatten
  def priors:IndexedSeq[Double]
  lazy val paramMap = models.map{_.paramMap}.foldLeft(mixtureParams){_++_}
  def mixtureParams=Map[(ParamName,Int),Any]()
  def updatedVec(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]):StdMixtureModel
  def updatedMat(p:ParamName,vec:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int]):StdMixtureModel
  def getOptParam(p:ParamName,index:Option[Int]):Option[IndexedSeq[Double]]
  def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]):StdMixtureModel
  def add(model:Model,prior:Double):StdMixtureModel
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
}
object Timing{
  def time[A](f: =>A):(A,Double)={
    val start = System.nanoTime : Double
    val ans = f
    val end = System.nanoTime : Double
    (ans,(end-start)/1000000.0)
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
  def zeroRate(piValues:IndexedSeq[Double],subModelIndex:Int)={
    new InvarLikelihoodModel(piValues,subModelIndex)
  }
}
/*
  Site class model with one n*m n*m matrix for n site classes of m-length alphabet.
  Subclasses could implement Covarion or THMM type models.
 */
object StdSiteClassModel{
  def apply(subModel:Seq[Model],priors:IndexedSeq[Double],modelIndex:Int,bigMat:IndexedSeq[IndexedSeq[Double]]) = {
    new StdSiteClassModel(subModel,priors,modelIndex,Some(bigMat),None)
  }
  def apply(subModel:Seq[Model],priors:IndexedSeq[Double],modelIndex:Int)={
    new StdSiteClassModel(subModel,priors,modelIndex,None,None)
  }
  def apply(models:Seq[Model]):StdSiteClassModel={
    new StdSiteClassModel(models,Vector.fill(models.length)(1.0/models.length),0,None,None)
  }
}

class ThmmSiteClassModel(realSwitchingModel:StdSiteClassModel,sigma:Double,alphabet:Alphabet,var bigMat:Option[IndexedSeq[IndexedSeq[Double]]],var myExp:Option[Exp]) extends StdModel{
  override def numClasses = realSwitchingModel.numClasses
  override lazy val paramMap=realSwitchingModel.paramMap.updated((Sigma,modelIndex),sigma)
  def modelIndex=realSwitchingModel.modelIndex
  lazy val mat = {
    bigMat getOrElse{
      val len = alphabet.matLength
      var ans = realSwitchingModel.mat
      for (i <- 0 until numClasses){
        for (j <- i+1 until numClasses){
          for (x <- 0 until len){
            val myI = i * len + x
            val myJ = j * len +x
            ans = ans.updated(myI,ans(myI).updated(myJ,sigma * pi()(myJ)))
            ans = ans.updated(myJ,ans(myJ).updated(myI,sigma * pi()(myI)))
          }
        }
      }
      ans = ans.fixDiag
      bigMat=Some(ans)
      ans
    }
  }
  def pi() = realSwitchingModel.pi
  override def pi(n:Node) = realSwitchingModel.pi(n)
  lazy val exp = myExp.getOrElse{val ans = DefaultExpFactory(mat); myExp = Some(ans); ans}
  def mySigma(p:ParamName,paramIndex:Option[Int])={ (p==Sigma && (paramIndex.isEmpty || paramIndex.get==modelIndex)) }
  def preserveCache(p:ParamName,paramIndex:Option[Int])={
    (p,paramIndex) match {
      case (Sigma,Some(modelIndex))=>true
      case (Sigma,None)=>true
      case _ => realSwitchingModel.preserveCache(p,paramIndex)
    }
  }
  def updatedVec(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int])={
    if (mySigma(p,paramIndex)){
      new ThmmSiteClassModel(realSwitchingModel,vec(0),alphabet,bigMat,myExp)
    }else {
      val cache = if (preserveCache(p,paramIndex)){ (bigMat,myExp)}else {(None,None)}
      new ThmmSiteClassModel(realSwitchingModel.updatedVec(p,vec,paramIndex),sigma,alphabet,cache._1,cache._2)
    }
  }
  def updatedMat(p:ParamName,vec:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int])= {
      val cache = if (preserveCache(p,paramIndex)){ (bigMat,myExp)}else {(None,None)}
      new ThmmSiteClassModel(realSwitchingModel.updatedMat(p,mat,paramIndex),sigma,alphabet,cache._1,cache._2)
  }
  def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int])={
    p match {
      case Sigma if (paramIndex.isEmpty || paramIndex.get==modelIndex) => new ThmmSiteClassModel(realSwitchingModel,vec(0),alphabet,None,None)
      case any => {
        val cache = if (preserveCache(p,paramIndex)){ (bigMat,myExp)}else {(None,None)}
        new ThmmSiteClassModel(realSwitchingModel.setOptParam(p,vec,paramIndex),sigma,alphabet,cache._1,cache._2)
      }
    }
  }
  def getOptParam(p:ParamName,paramIndex:Option[Int])={
     if (p==Sigma && (paramIndex.isEmpty || paramIndex.get==modelIndex)){
       Some(Vector(sigma))
     }else {
       realSwitchingModel.getOptParam(p,paramIndex)
     }
   }
 
  def updatedRate(r:Double)={
    new ThmmSiteClassModel(realSwitchingModel.updatedRate(r),sigma,alphabet,None,None)
  }

}

abstract class AbstractSiteClassModel extends StdMixtureModel with CombiningMatrixModel{
  def subModel:Seq[Model]
  override def models = subModel
  def priors:IndexedSeq[Double]
  def modelIndex:Int
  var cachedMat:Option[IndexedSeq[IndexedSeq[Double]]]
  var cachedExp:Option[Exp]
  override val numClasses = priors.length
  def mat = intMat
  private lazy val intMat = cachedMat getOrElse {
    val ans = combineMatrices(subModel.map{_.mat},pi)
    cachedMat = Some(ans)
    ans
  }
  lazy val pi = subModel.map{_.pi}.zip(priors).map{t=> t._1.map{_*t._2}}.flatten.toIndexedSeq
  lazy val exp = cachedExp.getOrElse{val ans = DefaultExpFactory(mat); cachedExp = Some(ans); ans}
}

trait CombiningMatrixModel{
  def combineMatrices(mats:Seq[Matrix],pi:IndexedSeq[Double])={
    mats.foldLeft[IndexedSeq[IndexedSeq[Double]]]( Vector[Vector[Double]]() ){(m,m2)=>m addClass m2}.normalise(pi)
  }
}
class StdSiteClassModel(val subModel:Seq[Model],val priors:IndexedSeq[Double],val modelIndex:Int=0,var cachedMat:Option[IndexedSeq[IndexedSeq[Double]]],var cachedExp:Option[Exp]) extends AbstractSiteClassModel{
  override val numClasses = subModel.foldLeft(0){_+_.numClasses}
  override lazy val paramMap = subModel.map{_.paramMap}.reduceLeft{_++_}.updated((MixturePrior,modelIndex),priors)
  def add(model:Model,prior:Double)={
    new StdSiteClassModel(subModel :+ model,priors.map{_*(1.0-prior)} :+ prior,modelIndex,None,None)
  }
  def preserveCache(p:ParamName,paramIndex:Option[Int])={
    (p,paramIndex) match {
      case (MixturePrior,Some(modelIndex))=>false
      case (MixturePrior,None)=>false
      case _ => false 
    }
  }
  def updatedVec(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]):StdSiteClassModel={
    p match {
      case MixturePrior if (paramIndex.isEmpty || paramIndex.get==modelIndex) => new StdSiteClassModel(subModel,vec,modelIndex,None,None)//cachedMat,cachedExp)
      case any => new StdSiteClassModel(subModel.map{_.updatedVec(p,vec,paramIndex)},priors,modelIndex,None,None)
    }
  }
  def updatedMat(p:ParamName,vec:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int]):StdSiteClassModel= {
     new StdSiteClassModel(subModel.map{_.updatedMat(p,vec,paramIndex)},priors,modelIndex,None,None)//,cachedMat,cachedExp)
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
  def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]):StdSiteClassModel={
    p match {
      case MixturePrior if (paramIndex.isEmpty || paramIndex.get==modelIndex) => new StdSiteClassModel(subModel,MixturePrior.getReal(vec),modelIndex,None,None)//cachedMat,cachedExp)
      case any => new StdSiteClassModel(subModel.map{_.setOptParam(p,vec,paramIndex)},priors,modelIndex,None,None)
    }
  }
  def updatedRate(r:Double):StdSiteClassModel={error ("Rate currently unimplemented for StdSiteClassModel"); this}
}

/*
  Uses an explicit zeroed rate matrix for inclusion in a THMM/Covarion model
*/
class InvarLikelihoodModel(piValues:IndexedSeq[Double],override val subModelIndex:Int=0) extends StdModel{
  lazy val mat = EnhancedIndexedMatrix.zero(piValues.length)
  lazy val exp = new Exp{
    val intExp = EnhancedLinearMatrix.eye(piValues.length)
    def exp(bl:Double)={
      intExp
    }  
  }
  def pi = piValues
  def updatedRate(r:Double)={this}
  def updatedVec(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int])={
    val mI = paramIndex getOrElse subModelIndex
    (p,mI) match {
      case (Pi,`subModelIndex`) => new InvarLikelihoodModel(vec,subModelIndex) 
    }
  }
  def updatedMat(p:ParamName,vec:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int])= {
    this
  }
  def getOptParam(p:ParamName,modelIndex:Option[Int])={
    val mI = modelIndex.getOrElse(subModelIndex)
    (p,mI) match {
      case (Pi,`subModelIndex`)=> Some(Pi.getOpt(piValues))
      case (_,_) => None
    }
  }
  def setOptParam(p:ParamName,vec:IndexedSeq[Double],modelIndex:Option[Int])={
    val mI = modelIndex.getOrElse(subModelIndex)
    (p,mI) match {
      case (Pi,`subModelIndex`)=> new InvarLikelihoodModel(Pi.getReal(vec),subModelIndex)
    }
  }

}
class BasicLikelihoodModel(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]], override val rate:Double, val exp:Exp,override val subModelIndex:Int=0) extends StdModel{
  val mat = s.sToQ(piValues,rate)
  def pi = piValues
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
  def apply(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]],alpha:Double,numCat:Int,myParamIndex:Int,basicModels1:Model):GammaModel={
    apply(piValues,s,alpha,numCat,myParamIndex,Gamma.getDist(alpha,numCat).map{t=>basicModels1.updatedRate(t)})  
  }
  def apply(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]],alpha:Double,numCat:Int,myParamIndex:Int,models:Seq[Model]):GammaModel={
    new GammaModel(piValues,s,alpha,numCat,myParamIndex,models)
  }
}
class GammaModel(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]],alpha:Double,numCat:Int,myParamIndex:Int,override val models:Seq[Model],val modelIndex:Int=0,var cachedMat:Option[Matrix]=None,var cachedExp:Option[Exp]=None) extends Model with CombiningMatrixModel{
  val priors = Vector.fill(numCat)(1.0/numCat)
  lazy val pi = { models.map{_.pi}.zip(priors).map{t=> t._1.map{_*t._2}}.flatten.toIndexedSeq }
  def pis = models.map{_.pi}
  def mats = models.map{_.mat}
  def updatedRate(rate:Double)=error("Unimplemented")
  def paramMap = models.head.paramMap.updated((Gamma,modelIndex),alpha)

  lazy val numClasses = models.length

  lazy val mat = cachedMat.getOrElse{
    cachedMat=Some(combineMatrices(models.map{_.mat},pi))
    cachedMat.get
  }
  lazy val exp = cachedExp.getOrElse{
    cachedExp=Some(DefaultExpFactory(mat))
    cachedExp.get
  }

  def add(m:Model,prior:Double)={
    StdSiteClassModel(List(this,m),Vector(1.0-prior,prior),modelIndex+1)
  }

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

