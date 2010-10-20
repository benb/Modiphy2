package modiphy.opt
import modiphy.tree._
import modiphy.calc._
import modiphy.alignment._
import modiphy.model._

sealed trait ParamName{
  def apply(i:Int)=(this,i)
  def lower(i:Int):Double
  def upper(i:Int):Double 
  def getReal(v:IndexedSeq[Double]):Any
  def from(o:Optimizable):Option[Any] = {o(this) match {
    case None => None
    case Some(ans) => Some(getReal(ans))
  }}
  def getRand(current:IndexedSeq[Double])={
    import scala.util.Random._
    current.zipWithIndex.map{ t=>
      val l = lower(t._2)
      val u = upper(t._2)
      nextDouble() * (u-l) - l
    }
  }
}

trait TypedParamName[A] extends ParamName{
  def getReal(v:IndexedSeq[Double]):A
  override def from(o:Optimizable):Option[A] = {o(this) match {
    case None => None
    case Some(ans) => Some(getReal(ans))
  }}
}
case class VectorParamWrapper(s:SingleParamName) extends VectorParamName{
  def lower(i:Int)=s.lower(i)
  def upper(i:Int)=s.upper(i)
}
trait SingleParamName extends TypedParamName[Double]{
  def getOpt(v:Double)=Vector(v)
  def getReal(v:IndexedSeq[Double])=v(0)
  def <<(v:Double)=(this,getOpt(v))
}
trait VectorParamName extends TypedParamName[IndexedSeq[Double]]{
  def getOpt(v:IndexedSeq[Double])=v
  def getReal(v:IndexedSeq[Double])=v
  def <<(v:IndexedSeq[Double])=(this,getOpt(v))
}
trait MatrixParamName extends TypedParamName[IndexedSeq[IndexedSeq[Double]]]{
  def getOpt(v:IndexedSeq[IndexedSeq[Double]]):IndexedSeq[Double]
  def getReal(v:IndexedSeq[Double]):IndexedSeq[IndexedSeq[Double]]
  def <<(v:IndexedSeq[IndexedSeq[Double]])=(this,getOpt(v))
}

case object Pi extends PiParamName
case object MixturePrior extends PiParamName 
case object Rate extends SingleParamName{
  def lower(i:Int)=0.0
  def upper(i:Int)=10000.0
}
case object Sigma extends SingleParamName{
  def lower(i:Int)=0.0
  def upper(i:Int)=100.0
}
trait PiParamName extends VectorParamName{
  override def getReal(d:IndexedSeq[Double])={
    val exponentiated =  d.map{i=>math.exp(i)}
    val total = exponentiated.reduceLeft{_+_} + math.exp(0.0D)
    val ans = new scala.collection.immutable.VectorBuilder[Double]
    ans += math.exp(0.0D)/total 
    ans ++= (0 until d.length).map{i=> exponentiated(i)/total}
    ans.result
  }
  override def getOpt(d:IndexedSeq[Double])={
    val head=d.head
    d.tail.map{i=>math.log(i/head)}.zipWithIndex.map{t=> val (i,c)=t
      if (i < lower(c)){lower(c)}else if(i > upper(c)){upper(c)} else {i}
    }
  }
  def lower(i:Int)= -10.0
  def upper(i:Int)=10.0
  override def getRand(current:IndexedSeq[Double])={
    import scala.util.Random._
    import scala.math._
    val ans = current.zipWithIndex.map{ t=>
      nextDouble()
    }.map{log}.map{_*(nextInt(2)*2-1)}
    println("Rand Pi " + ans)
    ans
  }
}
case object S extends MatrixParamName{
  def lower(i:Int)=0.0
  def upper(i:Int)=50.0
  def getOpt(d:Seq[Seq[Double]]):IndexedSeq[Double]=getOpt(d.map{_.toIndexedSeq}.toIndexedSeq)
  def getOpt(d:IndexedSeq[IndexedSeq[Double]]):IndexedSeq[Double]={
    var initial = d(1)(0) match {
      case 1.0=>d
      case norm => d.map{_.map{_/norm}}
    }
    val ans = for (i<-1 until d.length) yield {
     initial = initial.tail 
     initial.head.take(i)
    }
    ans.flatten.tail
  }
  def getReal(d:IndexedSeq[Double])={
    val len = d.length+1
    val matsize = 1 + (math.sqrt(8*len+1).toInt -1)/2 //from triangular root
    val initial = List(1.0).iterator ++ d.iterator
    for (i<-0 until matsize) yield {
      {initial.take(i).toList ++ List.fill(matsize-i)(0.0)}.toIndexedSeq
    }
  }
}

case object BranchLengths extends SingleParamName{
  def lower(i:Int)=0.0
  def upper(i:Int)=100.0
}

case object Gamma extends SingleParamName{
  lazy val gamma = new GammaDist
  def lower(i:Int)=0.01
  def upper(i:Int)=1000.0
  def getDist(d:Double,nC:Int)={
    gamma(d,nC)
  }
}

sealed abstract class ParamMatcher{
  def apply(i:Int):Boolean
}
sealed abstract class NumParamMatcher extends ParamMatcher{
  def params:Seq[Int]
}
case class MatchP(i:Int) extends NumParamMatcher{
  def apply(j:Int)=(i==j)
  def params=List(i)
}
case object MatchAll extends ParamMatcher{
  def apply(j:Int)=true
}
case class MatchSet(s:IndexedSeq[Int]) extends NumParamMatcher{
  def apply(j:Int)={
    s contains j
  }
  def params = s
}

object MatchSet{
  def apply(l:Int*):MatchSet={
    MatchSet(l.toIndexedSeq)
  }
}

class JointOptModel(val l:IndexedSeq[OptModel]) extends Optimizable{

  val numBL = l.map{_.tree.getBranchLengths.length}
  val blIndices=numBL.take(l.length-1).foldLeft(List[Int](0)){(l,d)=>(l.head+d)::l}
  println("DEBUG blIndices" + blIndices)
  assert(blIndices.length == l.length)
  var loc = 0
  val blMap:Map[Int,(Int,OptModel)] = blIndices.foldLeft((Map[Int,(Int,OptModel)](),0)){(t,index)=> val (m,last)=t
    val m2 = (last until index).foldLeft(m){(m2,i)=>
      m.updated(i,((i-last),l(loc)))
    }
    loc=loc+1
    (m2,index)
  }._1

  def calcSeq = l.map{_.calcSeq}.flatten
  def toXML = <JointModel>
  {l.zipWithIndex.map{t=> <model index={t._2.toString}>t._1.toXML</model>}}
 </JointModel>

  
  def logLikelihood = l.map{_.logLikelihood}.reduceLeft(_+_)

 def update(p:ParamName,d:Double,paramIndex:ParamMatcher){
    p match {
      case BranchLengths => 
        paramIndex match {
          case MatchAll => cantHandle(p,d,paramIndex)
          case MatchP(i) => val (i2,opt) = blMap(i); opt.update(p,d,MatchP(i2))
          case MatchSet(s) => 
            val s2 = s.map{blMap}
            l.foreach{opt=>
              val s3= s2.filter{t=> t._2==opt}.map{t=>t._1}
              if (s3.size>0){
                opt.update(p,d,MatchSet(s3))
              }
            }
        }
        case other => l.foreach{opt=> opt.update(p,d,paramIndex)}
    }
  }

  def update(p:ParamName,value:IndexedSeq[Double],paramIndex:ParamMatcher){
    p match {
      case BranchLengths => 
        paramIndex match {
          case MatchAll => {
            val vIter = value.iterator
            val newValues = numBL.map{i=> vIter.take(i).toIndexedSeq}
            l.zip(newValues).foreach{t=> t._1.update(p,t._2,paramIndex)}
          }
          case MatchP(i) => val (i2,opt) = blMap(i); opt.update(p,value,MatchP(i2))
          case MatchSet(i)=> 
            val s2 = i.map{blMap}
            l.foreach{opt=>
              val s3= s2.filter{t=> t._2==opt}.map{t=>t._1}
              if (s3.size>0){
                opt.update(p,value,MatchSet(s3))
              } 
            } 

        }
        case other => l.foreach{opt=> opt.update(other,value,paramIndex)}
    }
  }

  def update(p:ParamName,value:IndexedSeq[IndexedSeq[Double]],paramIndex:ParamMatcher=MatchAll)(implicit m:Manifest[IndexedSeq[IndexedSeq[Double]]]){
    p match {
      case BranchLengths=> cantHandle(p,value,paramIndex)
      case other => l.foreach{ opt => opt.update(other,value,paramIndex)}
    }
  }

  def setOptParam(p:ParamName,value:IndexedSeq[Double],paramIndex:ParamMatcher){
    p match {
      case BranchLengths => update(BranchLengths,value,paramIndex)
      case other => l.foreach{ opt => 
        opt.setOptParam(other,value,paramIndex)
      }
    }
  }

  def getOptParam(p:ParamName,paramIndex:ParamMatcher=MatchAll):Option[IndexedSeq[Double]]=(p,paramIndex) match {
      case (BranchLengths,MatchAll) => {
        Some(
          l.map(_.getOptParam(p,paramIndex)).filterNot{_==None}.map{_.get}.flatten
        )
      }
      case (BranchLengths,MatchP(i)) => 
        val (i2,opt) = blMap(i); opt.getOptParam(BranchLengths,MatchP(i2))
      case (BranchLengths,MatchSet(i)) => 
        val s2 = i.map{blMap}
        Some(
        l.map{opt=>
          val s3= s2.filter{t=> t._2==opt}.map{t=>t._1}
          if (s3.size>0){
            opt.getOptParam(p,MatchSet(s3))
          }else {
            None
          }
        }.filterNot{_==None}.map{_.get}.flatten
        )
      case _ => l.view.map{_.getOptParam(p,paramIndex)}.find(_!=None).getOrElse(None)
    }
def cantHandle(p:ParamName,a:Any,paramIndex:ParamMatcher){
    println("Can't handle combination " + p + " " + paramIndex + " " + a)
  }

  def apply(t:(ParamName,Int)):Option[IndexedSeq[Double]] = l.view.map{_(t)}.find{_!=None}.getOrElse{None}
  def apply(p:ParamName):Option[IndexedSeq[Double]] = l.view.map{_(p)}.find{_!=None}.getOrElse{None}


}
class OptModel(var calc:LikelihoodCalc,var tree:Tree,aln:Alignment) extends Optimizable{

  def calcSeq=List(calc)
  def m = calc.model
  val myParams:List[(ParamName,Int)] ={ m.numberedParams ++ tree.getBranchLengths.zipWithIndex.map{t=>(BranchLengths,t._2)}}.toList


   def update(p:ParamName,d:Double,paramIndex:ParamMatcher){
    p match {
      case BranchLengths => 
        paramIndex match {
          case MatchAll => cantHandle(p,d,paramIndex)
          case MatchP(i) => tree = tree.setBranchLength(i,d);updatedTree()
          case MatchSet(s) => {tree = s.foldLeft(tree){(t,i)=>t.setBranchLength(i,d)};updatedTree()}
        }
      case s:SingleParamName => calc = calc.updatedSingle(s,d,paramIndex)
    }
  }

  def update(p:ParamName,value:IndexedSeq[Double],paramIndex:ParamMatcher){
    p match {
      case BranchLengths => 
        paramIndex match {
          case MatchAll => tree = tree.setBranchLengths(value);updatedTree()
          case MatchSet(i)=> tree = i.zip(value).foldLeft(tree){(t,iv)=> t.setBranchLength(iv._1,iv._2)};updatedTree()
          case MatchP(i) => update(p,value(0),paramIndex)
        }
      case v:VectorParamName => calc = calc.updatedVec(v,value,paramIndex)
      case s:SingleParamName => calc = calc.updatedSingle(s,value(0),paramIndex)
    }
  }

  def update(p:ParamName,value:IndexedSeq[IndexedSeq[Double]],paramIndex:ParamMatcher=MatchAll)(implicit m:Manifest[IndexedSeq[IndexedSeq[Double]]]){
    p match {
      case BranchLengths=> cantHandle(p,value,paramIndex)
      case m:MatrixParamName => calc = calc.updatedMat(m,value,paramIndex)
    }
  }

  def setOptParam(p:ParamName,value:IndexedSeq[Double],paramIndex:ParamMatcher){
    p match {
      case BranchLengths => update(BranchLengths,value,paramIndex)
      case _ => calc = calc.setOptParam(p,value,paramIndex)
    }
  }

  def getOptParam(p:ParamName,paramIndex:ParamMatcher=MatchAll)={
    (p,paramIndex) match {
      case (BranchLengths,MatchAll) => {
        println("Getting branch lengths " + tree.getBranchLengths)
        Some(tree.getBranchLengths)
      }
      case (BranchLengths,MatchP(i)) => Some(Vector(tree.getBranchLengths(i)))
      case (BranchLengths,MatchSet(i)) => Some(i.map{j=>tree.getBranchLengths(j)})
      case _ => m.getOptParam(p,paramIndex)
    }
  }
  def logLikelihood=calc.logLikelihood
  

  def updatedTree(){
    calc = calc updated tree
  }
  def cantHandle(p:ParamName,a:Any,paramIndex:ParamMatcher){
    println("Can't handle combination " + p + " " + paramIndex + " " + a)
  }

  def apply(t:(ParamName,Int)):Option[IndexedSeq[Double]] = m.getOptParam(t._1,MatchP(t._2))
  def apply(p:ParamName):Option[IndexedSeq[Double]]=m.getOptParam(p,MatchAll)
  def toXML =  calc.toXML


}
trait Optimizable{

  def update(p:ParamName,d:Double){update(p,d,MatchAll)}
  def update(p:ParamName,d:IndexedSeq[Double]){update(p,d,MatchAll)}
  def update(p:ParamName,d:IndexedSeq[IndexedSeq[Double]])(implicit m:Manifest[IndexedSeq[IndexedSeq[Double]]]){update(p,d,MatchAll)}

  def update(t:(ParamName,Int), value:Double){ update(t._1,value,MatchP(t._2))}
  def update(t:(ParamName,Int), value:IndexedSeq[Double]){ update(t._1,value,MatchP(t._2))}
  def update(t:(ParamName,Int), value:IndexedSeq[IndexedSeq[Double]])(implicit m:Manifest[IndexedSeq[IndexedSeq[Double]]]){ update(t._1,value,MatchP(t._2))}
  def update(p:ParamName,d:IndexedSeq[Double],m:ParamMatcher):Unit
  def update(p:ParamName,d:Double,m:ParamMatcher):Unit
  def update(p:ParamName,d:IndexedSeq[IndexedSeq[Double]],m:ParamMatcher)(implicit man:Manifest[IndexedSeq[IndexedSeq[Double]]]):Unit


  def apply(t:(ParamName,Int)):Option[IndexedSeq[Double]]
  def apply(p:ParamName):Option[IndexedSeq[Double]]
  def calcSeq:Seq[LikelihoodCalc]
  def getOptParam(p:ParamName,paramIndex:ParamMatcher):Option[IndexedSeq[Double]]
  def getOptParams(optParams:Seq[(ParamName,ParamMatcher)]):Seq[IndexedSeq[Double]]={
    optParams.map{t=> getOptParam(t._1,t._2).getOrElse(error("Can't find param " + t))}
  }
  def logLikelihood:Double
  def setOptParam(p:ParamName,values:IndexedSeq[Double],paramIndex:ParamMatcher):Unit
  import dr.math.{UnivariateFunction,UnivariateMinimum,MultivariateFunction,ConjugateDirectionSearch}
  trait ModiphyMultivariateFunction extends MultivariateFunction{
    def getBestParam:Array[Double]
  }
  def getFunc(optParams:Seq[(ParamName,ParamMatcher)]):Either[UnivariateFunction,ModiphyMultivariateFunction]={
  val start = getOptParams(optParams)
    val lengths = start.map{_.length}
    val startArray = start.flatten.toArray
    val numArguments = startArray.length
    println("Start " + startArray.toList)
    println("NumArg " + numArguments)

    if (numArguments > 1){
      var bestlnL = -1E100
      var bestParam:Array[Double]=startArray
      var count=0
      val func = new ModiphyMultivariateFunction{
        val lowerBound = optParams.zip(lengths).map{t=> val ((pName,pIndex),len) = t; (0 until len).map{i=> pName.lower(i)}}.flatten
        val upperBound = optParams.zip(lengths).map{t=> val ((pName,pIndex),len) = t; (0 until len).map{i=> pName.upper(i)}}.flatten


        def getLowerBound(i:Int)=lowerBound(i)
        def getUpperBound(i:Int)=upperBound(i)
        val getNumArguments = numArguments
        def evaluate(params:Array[Double])={
          if (params.find(_.isNaN).isDefined){
            error("optimiser supplied NaN! " + params.mkString(" "))
          }
          var p2 = params.toList
          val splitP = lengths.map{l=> val ans = p2.take(l).toIndexedSeq; p2 = p2.drop(l); ans}
          optParams.zip(splitP).foreach{ t=> val ((pName,pIndex),values)=t
            setOptParam(pName,values,pIndex)
          }
          val ans = try {
            logLikelihood
          }catch {
            case e:Exception=>Double.NaN
          }
          if (ans > bestlnL){
            bestParam = params
            bestlnL = ans
          }
          println("OPT " + count + " " + optParams + " " + params.toList + " " + ans)
          count=count+1
          if (ans.isNaN){
            1E100
          }else {
            -ans
          }
        }
        def getBestParam = bestParam
      }
      Right(func)
    }else {
      println("Single Func")
      val func = new UnivariateFunction{
        val (paramName,paramIndex) = optParams.head
        val getLowerBound = paramName.lower(0)
        val getUpperBound = paramName.upper(0)
        def evaluate(p:Double)={
          setOptParam(paramName,Vector(p),paramIndex)
          val ans = logLikelihood
          println("OPT " + p + " " + ans)
          -ans
        }
      }
      Left(func)
    }
  }

  
  def optimiseAll(list:ParamName*){optimise(list.map{(_,MatchAll)})}
  def optimise(list:(ParamName,ParamMatcher)*){optimise(list)}
  def optimise(optParams:Seq[(ParamName,ParamMatcher)])(implicit p:Manifest[Seq[(ParamName,ParamMatcher)]]){
    val eitherFunc = getFunc(optParams)
    val start = getOptParams(optParams)
    val startArray = start.flatten.toArray
    eitherFunc match {
      case Left(func)=>
        println("Opt single " + startArray(0))
        val search = new UnivariateMinimum
        val finalP = search.optimize(startArray(0),func,1E-4)
        func evaluate finalP
      case Right(func) => 
        println("Opt multiple")
        val search = new ConjugateDirectionSearch
        val ans = safeOpt(search,func,startArray)
        if (ans.isDefined){
          println("Best " + func.getBestParam.toList + " " + func.evaluate(func.getBestParam))
        }else {
          error("Can't Optimise")
        }

    }

      
  }

  def randomise(list:(ParamName,ParamMatcher)*){
    randomise(list)
  }
  def randomise(s:Seq[(ParamName,ParamMatcher)])(implicit p:Manifest[Seq[(ParamName,ParamMatcher)]]){
    val params = getOptParams(s)
    s.zip(params).foreach{ t => val(pT,start)=t
      setOptParam(pT._1, pT._1.getRand(start), pT._2)
    }
  }
  import scala.annotation.tailrec

  @tailrec
  final def safeOpt(search:ConjugateDirectionSearch,func:ModiphyMultivariateFunction,startArray:Array[Double]):Option[Seq[Double]]={
    val ans = 
    try {
      search.optimize(func,startArray.clone,1E-4,1E-3)
      Some(func.getBestParam.toList)
    }catch {
      case e:Exception=>
      println(e)
      None
    }
    if (ans.isDefined){
      ans
    }else if (func.getBestParam!=startArray){
      println("Restarting with " + func.getBestParam.mkString(" ") + "(startarray = " + startArray.mkString(" ") + ")")
      safeOpt(search,func,func.getBestParam)
    }else {
      None
    }
  }

  def toXML:scala.xml.Elem

}


