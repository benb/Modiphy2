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
  type Parameters = Map[ParamName,IndexedSeq[Double]]
  type CalcCache = Map[Symbol,Any]
  implicit def WrapSingle(p:SingleParamName)=VectorParamWrapper(p)
}
import Types._

trait Model{
  def params:List[(ParamName,Int)]
  def updatedVec(p:VectorParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]):Model
  def updatedSingle(p:SingleParamName,vec:IndexedSeq[Double],paramIndex:Option[Int])={
    updatedVec(VectorParamWrapper(p),vec,paramIndex)
  }
  def updatedMat(p:MatrixParamName,vec:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int]):Model
  def getOptParam(p:ParamName,index:Option[Int]=None):Option[IndexedSeq[Double]]
  def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]):Model
  def updatedRate(r:Double):Model
  def scale(scale:Double)=updatedRate(rate*scale)
  def numClasses:Int
  def mats:Seq[Matrix]
  def pis:Seq[IndexedSeq[Double]]
  def pi:IndexedSeq[Double]
  def mat:Matrix
  def pi(n:Node):IndexedSeq[Double]=pi
  val exp:Exp
  def priors:IndexedSeq[Double]
  val rate=1.0
  def apply(e:Edge)=exp(e.dist)
  def models:Seq[Model]=List(this)
  def add(m:Model,prior:Double):Model
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
  def updatedVec(p:VectorParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]):StdModel
  def updatedMat(p:MatrixParamName,vec:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int]):StdModel
  def updatedRate(d:Double):StdModel
  def getOptParam(p:ParamName,index:Option[Int]):Option[IndexedSeq[Double]]
  def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]):StdModel
  def numClasses = 1
  def add(m:Model,p:Double)={
    new StdSiteClassModel(List(this,m),Vector(1.0-p,p),3,None,None)
  }
  /*
  def u:Matrix
  def d:Vector
  def v:Matrix
  */
}
trait StdMixtureModel extends Model{
  override def models:Seq[Model]
  def mats=models.map{_.mats}.flatten
  val rateCorrection = 1.0
  def pis = models.map{_.pis}.flatten
  def priors:IndexedSeq[Double]
  lazy val paramMap = models.map{_.params}.flatten.distinct.toList ++ mixtureParams.keys
  def mixtureParams=Map[(ParamName,Int),Any]()
  def updatedVec(p:VectorParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]):StdMixtureModel
  def updatedMat(p:MatrixParamName,vec:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int]):StdMixtureModel
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
  def rerate(newRate:Double):Exp = {
    val outer = this
    new Exp{
      def exp(bl:Double)=outer.exp(bl*newRate)
    }
  }
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
  override lazy val params= (Sigma,modelIndex) :: realSwitchingModel.params 
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
  def updatedVec(p:VectorParamName,vec:IndexedSeq[Double],paramIndex:Option[Int])={
    if (mySigma(p,paramIndex)){
      new ThmmSiteClassModel(realSwitchingModel,vec(0),alphabet,bigMat,myExp)
    }else {
      val cache = if (preserveCache(p,paramIndex)){ (bigMat,myExp)}else {(None,None)}
      new ThmmSiteClassModel(realSwitchingModel.updatedVec(p,vec,paramIndex),sigma,alphabet,cache._1,cache._2)
    }
  }
  def updatedMat(p:MatrixParamName,vec:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int])= {
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
    new ThmmSiteClassModel(realSwitchingModel.updatedRate(r),sigma,alphabet,bigMat,myExp)
  }

}

abstract class AbstractSiteClassModel extends StdMixtureModel with CombiningMatrixModel{
  def subModel:Seq[Model]
 def priors:IndexedSeq[Double]
  def modelIndex:Int
  var cachedMat:Option[IndexedSeq[IndexedSeq[Double]]]
  var cachedExp:Option[Exp]
  override val numClasses = priors.length
  def mat = intMat
  private lazy val intMat = cachedMat getOrElse {
    val ans = combineMatrices(subModel.map{_.mat},pi,rate)
    cachedMat = Some(ans)
    ans
  }
  lazy val pi = subModel.map{_.pi}.zip(priors).map{t=> t._1.map{_*t._2}}.flatten.toIndexedSeq
  lazy val exp = cachedExp.getOrElse{val ans = DefaultExpFactory(mat); cachedExp = Some(ans); ans}
}

trait CombiningMatrixModel{
  def combineMatrices(mats:Seq[Matrix],pi:IndexedSeq[Double],rate:Double)={
      mats.foldLeft[IndexedSeq[IndexedSeq[Double]]]( Vector[Vector[Double]]() ){(m,m2)=>m addClass m2}.normalise(pi,rate)
  }
}
class StdSiteClassModel(val subModel:Seq[Model],val priors:IndexedSeq[Double],val modelIndex:Int=0,var cachedMat:Option[IndexedSeq[IndexedSeq[Double]]],var cachedExp:Option[Exp],override val rate:Double=1.0) extends AbstractSiteClassModel{
  val params = (MixturePrior,modelIndex)::subModel.map{_.params}.flatten.distinct.toList
  override val numClasses = subModel.foldLeft(0){_+_.numClasses}
  override lazy val models = {
    //println("rate " + rate + " rateCorrection " + rateCorrection)
    subModel.map{r=>r.scale(rateCorrection)}
  }

  override def toString = "Mixture Model (Priors:" + priors + ") (\n" + models.map{_.toString.lines.map{x=> "  " + x}.mkString("\n")}.mkString("\n") + "\n)"
 
  override val rateCorrection = rate / (subModel.zip(priors).map{t=>t._1.rate * t._2}.reduceLeft{_+_})
  
  def updatedRate(r:Double):StdSiteClassModel=new StdSiteClassModel(subModel,priors,modelIndex,cachedMat,cachedExp,r)

  def add(model:Model,prior:Double)={
    new StdSiteClassModel(subModel :+ model,priors.map{_*(1.0-prior)} :+ prior,modelIndex,None,None,rate)
  }
  def preserveCache(p:ParamName,paramIndex:Option[Int])={
    (p,paramIndex) match {
      case (MixturePrior,Some(modelIndex))=>false
      case (MixturePrior,None)=>false
      case _ => false 
    }
  }
  def updatedVec(p:VectorParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]):StdSiteClassModel={
    p match {
      case MixturePrior if (paramIndex.isEmpty || paramIndex.get==modelIndex) => new StdSiteClassModel(subModel,vec,modelIndex,None,None,rate)//cachedMat,cachedExp)
      case any => new StdSiteClassModel(subModel.map{_.updatedVec(p,vec,paramIndex)},priors,modelIndex,None,None,rate)
    }
  }
  def updatedMat(p:MatrixParamName,vec:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int]):StdSiteClassModel= {
     new StdSiteClassModel(subModel.map{_.updatedMat(p,vec,paramIndex)},priors,modelIndex,None,None,rate)//,cachedMat,cachedExp)
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
      case MixturePrior if (paramIndex.isEmpty || paramIndex.get==modelIndex) => new StdSiteClassModel(subModel,MixturePrior.getReal(vec),modelIndex,None,None,rate)//cachedMat,cachedExp)
      case any => new StdSiteClassModel(subModel.map{_.setOptParam(p,vec,paramIndex)},priors,modelIndex,None,None,rate)
    }
  }
}


object BasicLikelihoodModel{
  def apply(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]]):BasicLikelihoodModel={
    apply(piValues,s,0)
  }
  def apply(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]],subModelIndex:Int):BasicLikelihoodModel={
    apply(piValues,s,1.0,subModelIndex)
  }
  def apply(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]],rate:Double,subModelIndex:Int):BasicLikelihoodModel={
    new BasicLikelihoodModel(piValues,s,rate,DefaultExpFactory(s,piValues,rate),subModelIndex)
  }
  def zeroRate(piValues:IndexedSeq[Double],subModelIndex:Int)={
    new InvarLikelihoodModel(piValues,subModelIndex)
  }
  val cleanParams=Set[ParamName](Pi,S)
}

trait UsefulModelUtil extends StdModel{
  def getCache[A](s:Symbol,calc:()=>A):A={
    cache.getOrElse(s,{
        val ans = calc()
        cache=cache.updated(s,ans)
        ans
      }) match {
        case ans:A=> ans
        case _ => error("Cache bug!")
    }
  }

  def cleanParams:Set[ParamName]
  def parameters:Parameters
  def modelIndex:Int
  var cache:CalcCache
  override def apply(e:Edge)=exp(e.dist)
  
 def update(newParameters:(ParamName,IndexedSeq[Double])*):UsefulModelUtil

  def updatedDouble(p:SingleParamName,vec:Double,paramIndex:Option[Int])={
    if (params.contains(p) && (paramIndex.isEmpty || paramIndex.get == modelIndex)){
      println("Updating " + p + " " + paramIndex)
      update((p,p.getOpt(vec)))
    }else {this}
  }
  def updatedVec(p:VectorParamName,vec:IndexedSeq[Double],paramIndex:Option[Int])={
    if (params.contains(p) && (paramIndex.isEmpty || paramIndex.get == modelIndex)){
      println("Updating " + p + " " + paramIndex)
      p match {
        case VectorParamWrapper(p2) => updatedDouble(p2,vec(0),paramIndex)
        case x=> update((p,p.getOpt(vec)))
      }
    }else {this}
  }
  def updatedMat(p:MatrixParamName,vec:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int])= {
    if (params.contains(p) && (paramIndex.isEmpty || paramIndex.get == modelIndex)){
      update((p,p.getOpt(vec)))
    }else {this}
  }

}
trait UsefulMixtureModel extends UsefulModelUtil{

 lazy val rateCorrection = rate / (subModels.zip(priors).map{t=>t._1.rate * t._2}.reduceLeft{_+_})
 def subModels:Seq[UsefulModelUtil]
 override lazy val models=subModels.map{_.scale(rateCorrection)}

 def factory(parameters:Parameters,subModel:Seq[Model],cache:CalcCache):UsefulMixtureModel
 override def update(newParameters:(ParamName,IndexedSeq[Double])*):UsefulMixtureModel={
    val newModels = subModels.map{_.update(newParameters:_*)}
    val changedModels = newModels.zip(subModels).foldLeft(false){(bool,t)=> bool || !(t._1 eq t._2)} // check for reference equality

    val unClean = changedModels || newParameters.foldLeft(false){(bool,t)=>
      bool || cleanParams.contains(t._1)
    } || true //FIXME
    val newP = newParameters.foldLeft(parameters){(m,t)=>
      if(m contains t._1){
        m updated (t._1,t._2)
      }else {
        m
      }
    }
    if (!(changedModels) && (newParameters==parameters)){
      this
    }else if (!unClean){
      factory(newP,newModels,cache)
    }else {
      factory(newP,newModels,Map[Symbol,Any]())
    }
  }

  override def toString:String={
    "---\n"
    super.toString + ":\n" + 
    parameters.map{t=>
      t._1 + " " + t._1.getReal(t._2) + "\n" 
    } + "SubModels (\n" + subModels.map{_.toString.lines.map{x=> "  " + x}.mkString("\n")}.mkString("\n") + "\n)"
    "---\n"
  }
 
  
  lazy val myParams = parameters.keys.map{(_,modelIndex)}.toList
  lazy val params = myParams ++ subModels.map{_.params}.flatten.distinct

  def getOptParam(p:ParamName,paramIndex:Option[Int])={
    if ((paramIndex.isEmpty || paramIndex.get == modelIndex) && parameters.contains(p)){
      Some(parameters(p))
    }else {
      models.foldLeft[Option[IndexedSeq[Double]]](None){(ans,mod)=>
        if (ans.isEmpty){
          mod.getOptParam(p,paramIndex)
        }else {
          ans
        }
      }
    }
  }
  def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int])={
    if (params.contains(p) && (paramIndex.isEmpty || paramIndex.get == modelIndex)){
      update((p,vec))
    }else {this}
  }

  def updatedRate(r:Double)={update((Rate,Vector(r)))} //FIXME RESCALE EXP

}
trait UsefulModel extends UsefulModelUtil{
 override def update(newParameters:(ParamName,IndexedSeq[Double])*):UsefulModel={
    val unClean = newParameters.foldLeft(false){(bool,t)=>
      bool || cleanParams.contains(t._1)
    } || true //FIXME
    val newP = newParameters.foldLeft(parameters){(m,t)=>m updated (t._1,t._2)}
    if (! unClean){
      factory(newP,cache)
    }else {
      factory(newP,Map[Symbol,Any]())
    }
  }
  override def toString:String={
    "---\n"
    super.toString + ":\n" + 
    parameters.map{t=>
      t._1 + " " + t._1.getReal(t._2) + "\n"
    } + 
    "---\n"
  }
  def factory(parameters:Parameters,cache:CalcCache):UsefulModel
 
  lazy val params = parameters.keys.map{(_,modelIndex)}.toList

  def updatedRate(r:Double)={update((Rate,Vector(r)))} //FIXME RESCALE EXP

  def getOptParam(p:ParamName,paramIndex:Option[Int])={
    if ((paramIndex.isEmpty || paramIndex.get == modelIndex) && parameters.contains(p)){
      Some(parameters(p))
    }else {None}
  }
  def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int])={
    if (params.contains(p) && (paramIndex.isEmpty || paramIndex.get == modelIndex)){
      update((p,vec))
    }else {this}
  }

}

class BasicLikelihoodModel(val parameters:Parameters,val modelIndex:Int,var cache:CalcCache=Map[Symbol,Any]()) extends UsefulModel{
  def factory(parameters:Parameters,cache:CalcCache)=new BasicLikelihoodModel(parameters,modelIndex,cache)
  def this(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]], rate:Double, exp:Exp,subModelIndex:Int) ={
    this(Map(Pi->Pi.getOpt(piValues),S->S.getOpt(s),Rate->Rate.getOpt(rate)),subModelIndex)
  }
  val cleanParams = BasicLikelihoodModel.cleanParams 

  lazy val pi = Pi.getReal(parameters(Pi))
  lazy val s = S.getReal(parameters(S))
  lazy val mat = getCache[Matrix]('Q, {() => val ans = s.sToQ(pi,rate); ans}) 
  lazy val exp = getCache[Exp]('Exp, {() => DefaultExpFactory(mat)}) 
  override val rate = parameters(Rate)(0)

}
/*
  Uses an explicit zeroed rate matrix for inclusion in a THMM/Covarion model
*/
class InvarLikelihoodModel(val parameters:Parameters,val modelIndex:Int,var cache:CalcCache=Map[Symbol,Any]()) extends UsefulModel{
  def this (piValues:IndexedSeq[Double],modelIndex:Int) = this(Map(Pi->Pi.getOpt(piValues)),modelIndex,Map[Symbol,Any]())
  def factory(p:Parameters,cache:CalcCache)=new InvarLikelihoodModel(p,modelIndex,cache)
  val cleanParams = Set[ParamName]()
  lazy val pi = Pi.getReal(parameters(Pi))
  lazy val mat = getCache[Matrix]('Q, {() => EnhancedIndexedMatrix.zero(pi.length)})
  override val rate = 0.0
  lazy val exp = getCache[Exp]('Exp,{ () => 
    new Exp{
      val intExp = EnhancedLinearMatrix.eye(parameters(Pi).length + 1)
      def exp(bl:Double)={
        intExp
      }  
    }}
  )
}




object GammaModel{
  def apply(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]],alpha:Double,numCat:Int,rate:Double=1.0):GammaModel={
    apply(piValues,s,alpha,numCat,0,rate)
  }
  def apply(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]],alpha:Double,numCat:Int,myParamIndex:Int,rate:Double):GammaModel={
    apply(piValues,s,alpha,numCat,myParamIndex,BasicLikelihoodModel(piValues,s),rate)
  }
  def apply(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]],alpha:Double,numCat:Int,myParamIndex:Int,basicModels1:Model,rate:Double):GammaModel={
    apply(piValues,s,alpha,numCat,myParamIndex,Gamma.getDist(alpha,numCat).map{t=>basicModels1.updatedRate(t*rate)},rate)  
  }
  def apply(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]],alpha:Double,numCat:Int,myParamIndex:Int,models:Seq[Model],rate:Double):GammaModel={
    new GammaModel(piValues,s,alpha,numCat,myParamIndex,models,None,None,rate)
  }
}
class GammaModel(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]],alpha:Double,numCat:Int,modelIndex:Int,override val models:Seq[Model],var cachedMat:Option[Matrix]=None,var cachedExp:Option[Exp]=None,override val rate:Double=1.0) extends Model with CombiningMatrixModel{

  val params = Pi::S::Gamma::Nil map {(_,modelIndex)}

  val priors = Vector.fill(numCat)(1.0/numCat)
  lazy val pi = { models.map{_.pi}.zip(priors).map{t=> t._1.map{_*t._2}}.flatten.toIndexedSeq }
  def pis = models.map{_.pi}
  def mats = models.map{_.mat}
  def updatedRate(rate:Double)=GammaModel(piValues,s,alpha,numCat,modelIndex,models.head,rate) //TODO optimise this

  lazy val numClasses = models.length

  lazy val mat = cachedMat.getOrElse{
    cachedMat=Some(combineMatrices(models.map{_.mat},pi,rate))
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
  def updatedVec(p:VectorParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]) = {
    val mI = paramIndex.getOrElse(modelIndex)
    (p,mI) match {
      case (VectorParamWrapper(Gamma),`modelIndex`)=>{
        GammaModel(piValues,s,vec(0),numCat,modelIndex,models.head,rate)
      }
      case (Pi,`modelIndex`)=>{
        GammaModel(vec,s,alpha,numCat,modelIndex,rate)
      }
      case _ => this
    }
  }
  def updatedMat(p:MatrixParamName,vec:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int]) = {
    val mI = paramIndex.getOrElse(modelIndex)
    (p,mI) match {
      case (S,`modelIndex`)=>GammaModel(piValues,vec,alpha,numCat,modelIndex,rate)
      case _ => this
    }
  }

  def getOptParam(p:ParamName,paramIndex:Option[Int])={
    val mI = paramIndex.getOrElse(modelIndex)
    (p,mI) match {
      case (Pi,`modelIndex`)=> Some(Pi.getOpt(piValues))
      case (S,`modelIndex`)=> Some(S.getOpt(s))
      case (Gamma,`modelIndex`) => Some(Vector(alpha))
      case (_,_) => None
    }
  }
  def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int])={
    val mI = paramIndex.getOrElse(modelIndex)
    (p,mI) match {
      case (Pi,`modelIndex`)=> updatedVec(Pi,Pi.getReal(vec),paramIndex)
      case (S,`modelIndex`)=> updatedMat(S,S.getReal(vec),paramIndex)
      case (Gamma,`modelIndex`) => {
        updatedVec(Gamma,vec,paramIndex)
      }
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


