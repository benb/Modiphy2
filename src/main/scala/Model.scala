package modiphy.model
import modiphy.math.EnhancedMatrix._
import modiphy.tree._
import modiphy.calc._
import modiphy.util._
import scala.collection.LinearSeq
import modiphy.opt._
import modiphy.alignment.{Alphabet,Alignment}
import modiphy.math._

object Types{
  type Matrix=IndexedSeq[IndexedSeq[Double]]
  type LinearMatrix = scala.collection.immutable.LinearSeq[scala.collection.immutable.LinearSeq[Double]]
  type Parameters = Map[ParamName,IndexedSeq[Double]]
  type CalcCache = Map[Symbol,Any]
  val cleanCache = Map[Symbol,Any]()
  implicit def WrapParam(p:Parameters)=new ParametersWrapper(p)
  implicit def MakeParamMap(l:List[(ParamName,IndexedSeq[Double])])=l.toMap
}
import Types._
class ParametersWrapper(par:Parameters){
  def viewParam(p:SingleParamName)=p.getReal(par(p))
  def viewParam(p:VectorParamName)=p.getReal(par(p))
  def viewParam(p:MatrixParamName)=p.getReal(par(p))
  def updatedParam(p:SingleParamName,d:Double)=par.updated(p,p.getOpt(d))
  def updatedParam(p:VectorParamName,v:IndexedSeq[Double])=par.updated(p,p.getOpt(v))
  def updatedParam(p:MatrixParamName,m:IndexedSeq[IndexedSeq[Double]])=par.updated(p,p.getOpt(m))
}

/**
 Base Model type. For flexibility all Models can be accessed as either a
 "Mixture" model (as a series of rate matrices + probabilities) or a single
 large rate matrix. 
*/
trait Model{
  def likelihoodCalc(t:Tree,aln:Alignment):LikelihoodCalc[Model]
  def applies(p:ParamName,pM:ParamMatcher):Boolean
  def specialUpdate:PartialFunction[(ParamName,IndexedSeq[Double]),Option[Model]]={case _ => None}
  def params:scala.collection.Set[ParamName]
  def numberedParams:scala.collection.Set[(ParamName,Int)]
  def updatedVec(p:VectorParamName,vec:IndexedSeq[Double],paramIndex:ParamMatcher):Model
  def updatedSingle(p:SingleParamName,d:Double,paramIndex:ParamMatcher):Model
  def updatedMat(p:MatrixParamName,vec:IndexedSeq[IndexedSeq[Double]],paramIndex:ParamMatcher):Model
  def getOptParam(p:ParamName,index:ParamMatcher=MatchAll):Option[IndexedSeq[Double]]
  def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:ParamMatcher):Model
  def updatedRate(r:Double):Model
  def scale(scale:Double)=updatedRate(rate*scale)
  def numClasses:Int
  def mats:Seq[Matrix]
  def pis:Seq[IndexedSeq[Double]]
  def pi:IndexedSeq[Double]
  def mat:Matrix
  val exp:Exp
  def priors:IndexedSeq[Double]
  def rate:Double
  def apply(tree:Tree,node:NonRoot)=exp(tree.branchLength(node))
  def models:Seq[Model]=List(this)
  def add(m:Model,prior:Double,modelIndex:Int):Model
  def update(newParameters:(ParamName,ParamMatcher,IndexedSeq[Double])*):Model
  protected def flatPriors(x:Int) = Vector.fill(x){1.0/x}
  protected def repeatedPi(pi:IndexedSeq[Double],n:Int)={
    val inPi = pi.map{_/n}
    (0 until n).map{i=>inPi}.flatten
  }
}

trait StdModel extends Model{
  val subModelIndex=0
  def numClasses = 1
  def add(m:Model,p:Double,modelIndex:Int)={ new StdSiteClassModel((MixturePrior << Vector(1.0-p,p))::Nil,modelIndex,List(this,m)) }
  def priors=Vector(1.0)
  def mats=List(mat)
  def pis =List(pi)
  /*
  def u:Matrix
  def d:Vector
  def v:Matrix
  */
}
trait StdMixtureModel extends Model{
  def rateCorrection:Double
}


/*
  Site class model with one n*m n*m matrix for n site classes of m-length alphabet.
  Subclasses could implement Covarion or THMM type models.
 */
object StdSiteClassModel{
  def apply(subModel:Seq[Model],priors:IndexedSeq[Double],modelIndex:Int,bigMat:IndexedSeq[IndexedSeq[Double]]) = {
    new StdSiteClassModel((MixturePrior << priors)::(Rate << 1.0)::Nil,modelIndex,subModel)
  }
  def apply(subModel:Seq[Model],priors:IndexedSeq[Double],modelIndex:Int)={
    new StdSiteClassModel((MixturePrior << priors)::(Rate << 1.0)::Nil,modelIndex,subModel)
  }
  def apply(models:Seq[Model]):StdSiteClassModel={
    val priors = Vector.fill(models.length)(1.0/models.length)
    new StdSiteClassModel((MixturePrior << priors)::(Rate << 1.0)::Nil,0,models)
  }
}

class ThmmSiteClassModel(val parameters:Parameters,val modelIndex:Int,realSwitchingModel:Model,len:Int,var cache:CalcCache=cleanCache) extends UsefulWrappedModel with UsefulSingleModel{
  override val numClasses = realSwitchingModel.numClasses
  lazy val mat = getCache[Matrix]('Q,{()=>
      val sigma = parameters viewParam Sigma
      var ans = realSwitchingModel.mat
      for (i <- 0 until numClasses){
        for (j <- i+1 until numClasses){
          for (x <- 0 until len){
            val myI = i * len + x
            val myJ = j * len +x
            ans = ans.updated(myI,ans(myI).updated(myJ,sigma * pi(myJ)))
            ans = ans.updated(myJ,ans(myJ).updated(myI,sigma * pi(myI)))
          }
        }
      }
      ans = ans.fixDiag
      ans
  })
  def wrapped = List(realSwitchingModel)
  def cleanParams = params
  def pi = realSwitchingModel.pi
  override def updatedRate(r:Double)=factory(parameters.updatedParam(Rate,r),wrapped,cache) 
  def factory(parameters:Parameters,realSwitchingModel:Seq[Model],cache:CalcCache)=new ThmmSiteClassModel(parameters,modelIndex,realSwitchingModel.head,len,cache)
}


class StdSiteClassModel(val parameters:Parameters,val modelIndex:Int,val subModels:Seq[Model],var cache:CalcCache=cleanCache) extends UsefulMixtureModel{
//  def this(val subModel:Seq[Model],val priors:IndexedSeq[Double],val modelIndex:Int=0,
 //   var cachedMat:Option[IndexedSeq[IndexedSeq[Double]]],var cachedExp:Option[Exp],override val rate:Double=1.0)
 //  = 
 lazy val pi = subModels.zip(priors).map{t=>t._1.pi.map{_*t._2}}.flatten.toIndexedSeq
 lazy val priors = parameters viewParam MixturePrior
 val cleanParams=params
 def factory(parameters:Parameters,subModels:Seq[Model],cache:CalcCache)=new StdSiteClassModel(parameters,modelIndex,subModels,cache)//TODO optimise
}

class GammaModel(val parameters:Parameters,val modelIndex:Int, var cache:CalcCache=cleanCache,numCat:Int) extends UsefulMixtureModel{
  override lazy val numClasses = numCat
  def priors = Vector.fill(numClasses){1.0/numClasses}
  lazy val pi = {
    val smallPi = parameters viewParam Pi map {_/numClasses}
    (0 until numClasses).map{i=> smallPi}.flatten[Double]
  }
  override def add(model:Model,prior:Double,newModelIndex:Int)={
    val priors = Vector(1.0-prior,prior)
    new StdSiteClassModel((MixturePrior << priors) :: (Rate << 1.0) :: Nil,newModelIndex,List(this,model))
  }
  lazy val subS = parameters viewParam S
  override lazy val subModels = getCache[Seq[Model]]('Models,{()=>
    base.exp.instantiated // hack to make sure Exp is constructed and so can be cheaply copied
    Gamma.getDist(parameters viewParam Gamma,numClasses).map{ subrate => 
      base.updatedSingle(Rate,subrate * rate,MatchAll)
    }
  })
  lazy val base = getCache[Model]('Base,{()=>{
    new BasicLikelihoodModel(parameters.filter{t=> t._1==Pi || t._1==S || t._1==Rate}, modelIndex,cleanCache)
  }})
  val cleanParams=params
  def factory(parameters:Parameters,subModels:Seq[Model],cache:CalcCache)=new GammaModel(parameters,modelIndex,cache,numClasses)//TODO optimise

  override def specialUpdate:PartialFunction[(ParamName,IndexedSeq[Double]),Option[Model]] = {
    ({
      case (Gamma,vec)=>Some(new GammaModel(parameters updated (Gamma,vec),modelIndex,cache-'Models,numCat))
    }:PartialFunction[(ParamName,IndexedSeq[Double]),Option[Model]]).orElse(super.specialUpdate)
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

trait UsefulModelUtil extends Model{

  def applies(p:ParamName,paramIndex:ParamMatcher):Boolean 
  def appliesToMe(p:ParamName,paramIndex:ParamMatcher):Boolean = {
    parameters.contains(p) && paramIndex(modelIndex)
  }
  def cached(s:Symbol)=cache contains s
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
  override def specialUpdate:PartialFunction[(ParamName,IndexedSeq[Double]),Option[Model]] = {
    case _ => None
  }

  lazy val exp = getCache[Exp]('Exp,{()=>DefaultExpFactory(mat)})
  def params = parameters.keySet.toSet
  def numberedParams = params map {p => (p,modelIndex)}

  def cleanParams:Set[ParamName]
  def parameters:Parameters
  def modelIndex:Int
  def rate = parameters viewParam Rate
  var cache:CalcCache

  
  def update(newParameters:(ParamName,ParamMatcher,IndexedSeq[Double])*):Model

  def updatedSingle(p:SingleParamName,vec:Double,paramIndex:ParamMatcher)={
    if (applies(p,paramIndex)){
      update((p,paramIndex,p.getOpt(vec)))
    }else {this}
  }
  def updatedVec(p:VectorParamName,vec:IndexedSeq[Double],paramIndex:ParamMatcher)={
    if (applies(p,paramIndex)){
      p match {
        case x=> update((p,paramIndex,p.getOpt(vec)))
      }
    }else {this}
  }
  def updatedMat(p:MatrixParamName,vec:IndexedSeq[IndexedSeq[Double]],paramIndex:ParamMatcher)= {
    if (applies(p,paramIndex)){
      update((p,paramIndex,p.getOpt(vec)))
    }else {this}
  }

  def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:ParamMatcher)={
    if (applies(p,paramIndex)){
      update((p,paramIndex,vec))
    }else {this}
  }
}

trait UsefulWrappedModel extends UsefulModelUtil{
  def wrapped:Seq[Model]
  override def applies(p:ParamName,pM:ParamMatcher)={
    appliesToMe(p,pM) || wrapped.foldLeft(false){(bool,model)=> bool || model.applies(p,pM)}
  }
  def factory(parameters:Parameters,models:Seq[Model],cache:CalcCache):UsefulWrappedModel
  override def update(myNewP:(ParamName,ParamMatcher,IndexedSeq[Double])*):Model={
    val newParameters = myNewP.filter{t=>appliesToMe(t._1,t._2)}.map{t=>(t._1,t._3)}
    val special = newParameters.foldLeft[Option[Model]](Some(this)){(optM,p)=>
      if (optM.isEmpty){None}else {optM.get.specialUpdate(p)}     
    }
    special.getOrElse{
      val newModels = wrapped.map{_.update(myNewP:_*)}

      val changedModels = newModels.zip(wrapped).foldLeft(false){(bool,t)=> bool || !(t._1 eq t._2)} // check for reference equality

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
    factory(newP,newModels,cleanCache)
  }
}
  }

  override lazy val params = wrapped.map{_.params}.foldLeft(super.params){(m,sub)=>m ++ sub}
  def add(model:Model,prior:Double,newModelIndex:Int = modelIndex):Model={
    val priors = Vector(1.0-prior,prior)
    new StdSiteClassModel((MixturePrior << priors) :: (Rate << 1.0) :: Nil,newModelIndex,List(this,model))
  }
  
  def getOptParam(p:ParamName,paramIndex:ParamMatcher)={
    if (!applies(p,paramIndex)){
        None
    }else if (appliesToMe(p,paramIndex)){
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


}


trait UsefulMixtureModel extends UsefulWrappedModel{

 def likelihoodCalc(t:Tree,aln:Alignment) = new MixtureLikelihoodCalc(t,aln,this)
  override def add(model:Model,prior:Double,newModelIndex:Int = modelIndex)={
    if (modelIndex==newModelIndex){
      val newPriors = priors.map{_*(1.0-prior)} :+ prior
      val newParam = parameters updatedParam (MixturePrior,newPriors)
      factory(newParam,subModels :+ model,cleanCache)
    }else {
      super.add(model,prior,newModelIndex)
    }
  }

 lazy val mat = getCache[Matrix]('Q,{()=>combineMatrices(subModels.map{_.mat},pi,rate)})
  def combineMatrices(mats:Seq[Matrix],pi:IndexedSeq[Double],rate:Double)={
    mats.foldLeft[IndexedSeq[IndexedSeq[Double]]]( Vector[Vector[Double]]() ){(m,m2)=>m addClass m2}.normalise(pi,rate)
  }
  override lazy val numClasses = subModels.foldLeft(0){_+_.numClasses}
  lazy val rateCorrection = rate / (subModels.zip(priors).map{t=>t._1.rate * t._2}.reduceLeft{_+_})
  def subModels:Seq[Model]
  def wrapped = subModels
  override lazy val models=subModels.map{
    _.scale(rateCorrection)
  }
  def pis = models.map{_.pi}
  def mats = models.map{_.mat}

  def factory(parameters:Parameters,subModel:Seq[Model],cache:CalcCache):UsefulMixtureModel
  override def toString:String={
    "---\n"
    super.toString + ":\n" + 
    parameters.map{t=>
      t._1 + " " + t._1.getReal(t._2) + "\n" 
    } + "SubModels (\n" + subModels.map{_.toString.lines.map{x=> "  " + x}.mkString("\n")}.mkString("\n") + "\n)" + 
    "---\n"
  }
 
  

    def updatedRate(r:Double)={update((Rate,MatchP(modelIndex),Vector(r)))} //FIXME RESCALE EXP

}
trait UsefulSingleModel extends UsefulModelUtil{
 def priors = Vector(1.0)
 def pis = List(pi)
 def numClasses =1
 def mats = List(mat)
 def likelihoodCalc(t:Tree,aln:Alignment) = new SimpleLikelihoodCalc(t,aln,this)
}
trait UsefulModel extends UsefulSingleModel{
 def add(m:Model,prior:Double,addedIndex:Int) = {
   val priors = Vector(1.0-prior,prior)
   new StdSiteClassModel((MixturePrior << priors)::(Rate << 1.0)::Nil,addedIndex,List(this,m))
 }

 def applies(p:ParamName,paramIndex:ParamMatcher)=appliesToMe(p,paramIndex)

  lazy val mat = getCache[Matrix]('Q, {() => val ans = (parameters viewParam S).sToQ(pi,rate); ans}) 
 
  override def update(newP:(ParamName,ParamMatcher,IndexedSeq[Double])*):Model={
    val newParameters = newP.filter(t=>appliesToMe(t._1,t._2)).map{t=>(t._1,t._3)}
    val special = newParameters.foldLeft[Option[Model]](Some(this)){(optM,p)=>
      if (optM.isEmpty){None}else {optM.get.specialUpdate(p)}     
    }
    special.getOrElse{
      val unClean = newParameters.foldLeft(false){(bool,t)=>
        bool || cleanParams.contains(t._1)
      } || true //FIXME
      val newP = newParameters.foldLeft(parameters){(m,t)=> if (m contains t._1){m updated (t._1,t._2)}else {m}}
      if (newP==parameters){
        this
      }
      else if (! unClean){
        factory(newP,cache)
      }else {
        factory(newP,cleanCache)
      }
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
 
 override def updatedRate(r:Double)={
   exp.instantiated
    val newCache = if(cached('Exp)){
      Map[Symbol,Any]('Exp -> (exp))
    }else {
      cleanCache
    }
    factory(parameters.updated(Rate,Rate.getOpt(r)),newCache)
 }

 

  def getOptParam(p:ParamName,paramIndex:ParamMatcher)={
    if (applies(p,paramIndex)){
      Some(parameters(p))
    }else {None}
  }

}

class BasicLikelihoodModel(val parameters:Parameters,val modelIndex:Int,var cache:CalcCache=cleanCache) extends UsefulModel{
  def factory(parameters:Parameters,cache:CalcCache)=new BasicLikelihoodModel(parameters,modelIndex,cache)
  def this(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]], rate:Double, exp:Exp,subModelIndex:Int) ={
    this(Map(Pi->Pi.getOpt(piValues),S->S.getOpt(s),Rate->Rate.getOpt(rate)),subModelIndex)
  }
  val cleanParams = BasicLikelihoodModel.cleanParams 

  lazy val pi = Pi.getReal(parameters(Pi))
  lazy val s = S.getReal(parameters(S))
  override val rate = parameters(Rate)(0)

}
/*
  Uses an explicit zeroed rate matrix for inclusion in a THMM/Covarion model
*/
class InvarLikelihoodModel(val parameters:Parameters,val modelIndex:Int,var cache:CalcCache=cleanCache) extends UsefulModel{
  def this (piValues:IndexedSeq[Double],modelIndex:Int) = this(Map(Pi->Pi.getOpt(piValues)),modelIndex,cleanCache)
  def factory(p:Parameters,cache:CalcCache)=new InvarLikelihoodModel(p,modelIndex,cache)
  val cleanParams = Set[ParamName]()
  lazy val pi = Pi.getReal(parameters(Pi))
  override lazy val mat = getCache[Matrix]('Q, {() => EnhancedIndexedMatrix.zero(pi.length)})
  override def rate = 0.0
  override lazy val exp = getCache[Exp]('Exp,{ () => 
    new Exp{
      val intExp = EnhancedLinearMatrix.eye(parameters(Pi).length + 1)
      def exp(bl:Double)={
        intExp
      }  
    }}
  )
}




object GammaModel{
  def apply(piValues:IndexedSeq[Double],s:IndexedSeq[IndexedSeq[Double]],alpha:Double,numCat:Int,myParamIndex:Int=0,rate:Double=1.0):GammaModel={
    new GammaModel(
      (Pi << piValues) :: (S << s) :: (Gamma << alpha) :: (Rate << rate) :: Nil,
      myParamIndex,
      cleanCache,
      numCat
    )
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


