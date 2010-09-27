package modiphy.opt
import modiphy.tree._
import modiphy.calc._
import modiphy.alignment._
import modiphy.model._

trait ParamName{
  def apply(i:Int)=(this,i)
  def lower(i:Int):Double
  def upper(i:Int):Double 
  def getReal(v:IndexedSeq[Double]):Any
}
case class VectorParamWrapper(s:SingleParamName) extends VectorParamName{
  def lower(i:Int)=s.lower(i)
  def upper(i:Int)=s.upper(i)
}
trait SingleParamName extends ParamName{
  def getOpt(v:Double)=Vector(v)
  def getReal(v:IndexedSeq[Double])=v(0)
  def <<(v:Double)=(this,getOpt(v))
}
trait VectorParamName extends ParamName{
  def getOpt(v:IndexedSeq[Double])=v
  def getReal(v:IndexedSeq[Double])=v
  def <<(v:IndexedSeq[Double])=(this,getOpt(v))
}
trait MatrixParamName extends ParamName{
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
    d.tail.map{i=>math.log(i/head)}
  }
  def lower(i:Int)= -10.0
  def upper(i:Int)=10.0
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


class OptModel[A <: Model](var calc:LikelihoodCalc[A],var tree:Tree,aln:Alignment){
  def m = calc.model
  val myParams:List[(ParamName,Int)] ={ m.numberedParams ++ tree.getBranchLengths.zipWithIndex.map{t=>(BranchLengths,t._2)}}.toList


  def update(p:ParamName,d:Double){update(p,d,MatchAll)}
  def update(p:ParamName,d:IndexedSeq[Double]){update(p,d,MatchAll)}
  def update(p:ParamName,d:IndexedSeq[IndexedSeq[Double]])(implicit m:Manifest[IndexedSeq[IndexedSeq[Double]]]){update(p,d,MatchAll)}
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

  def update(t:(ParamName,Int), value:Double){ update(t._1,value,MatchP(t._2))}
  def update(t:(ParamName,Int), value:IndexedSeq[Double]){ update(t._1,value,MatchP(t._2))}
  def update(t:(ParamName,Int), value:IndexedSeq[IndexedSeq[Double]])(implicit m:Manifest[IndexedSeq[IndexedSeq[Double]]]){ update(t._1,value,MatchP(t._2))}

  def apply(t:(ParamName,Int)):Option[IndexedSeq[Double]] = m.getOptParam(t._1,MatchP(t._2))
  def apply(p:ParamName):Option[IndexedSeq[Double]]=m.getOptParam(p,MatchAll)

  def optimiseAll(list:ParamName*){optimiseSeq(list.map{(_,MatchAll)})}
  def optimise(list:(ParamName,ParamMatcher)*){optimiseSeq(list)}
  def optimiseSeq(optParams:Seq[(ParamName,ParamMatcher)]){
    val start = optParams.map{t=> getOptParam(t._1,t._2).getOrElse(error("Can't find param " + t))}
    val lengths = start.map{_.length}
    val startArray = start.flatten.toArray
    val numArguments = startArray.length
    println("Start " + startArray.toList)

    import dr.math.{UnivariateFunction,UnivariateMinimum,MultivariateFunction,ConjugateDirectionSearch}
    if (numArguments > 1){
      var bestlnL = -1E100
      var bestParam:Array[Double]=startArray
      val func = new MultivariateFunction{
        val lowerBound = optParams.zip(lengths).map{t=> val ((pName,pIndex),len) = t; (0 until len).map{i=> pName.lower(i)}}.flatten
        val upperBound = optParams.zip(lengths).map{t=> val ((pName,pIndex),len) = t; (0 until len).map{i=> pName.upper(i)}}.flatten


        def getLowerBound(i:Int)=lowerBound(i)
        def getUpperBound(i:Int)=upperBound(i)
        val getNumArguments = numArguments
        def evaluate(params:Array[Double])={
          var p2 = params.toList
          val splitP = lengths.map{l=> val ans = p2.take(l).toIndexedSeq; p2 = p2.drop(l); ans}
          optParams.zip(splitP).foreach{ t=> val ((pName,pIndex),values)=t
            setOptParam(pName,values,pIndex)
          }
          val ans = logLikelihood
          if (ans > bestlnL){
            bestParam = params
            bestlnL = ans
          }
          println("OPT " + optParams + " " + params.toList + " " + ans)
          if (ans.isNaN){
            1E100
          }else {
            -ans
          }
        }
      }
      val search = new ConjugateDirectionSearch
      search.optimize(func,startArray,1E-4,1E-3)
      println("Best " + bestParam.toList + " " + func.evaluate(bestParam))
    }else {
      val func = new UnivariateFunction{
        val (paramName,paramIndex) = optParams.head
        val getLowerBound = paramName.lower(0)
        val getUpperBound = paramName.upper(0)
        def evaluate(p:Double)={
          setOptParam(paramName,Vector(p),paramIndex)
          val ans = logLikelihood
          println(p + " " + ans)
          -ans
        }
      }
      val search = new UnivariateMinimum
      val finalP = search.optimize(startArray(0),func,1E-4)
      func evaluate finalP
    }
  }
}





