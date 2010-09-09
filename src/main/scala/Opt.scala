package modiphy.opt
import modiphy.tree._
import modiphy.alignment._
import modiphy.model._

trait ParamName{
  def apply(i:Int)=(this,i)
  def lower(i:Int):Double
  def upper(i:Int):Double 
}

case object Pi extends PiParamName
case object MixturePrior extends PiParamName 
trait PiParamName extends ParamName{
  def getReal(d:IndexedSeq[Double])={
    val exponentiated =  d.map{i=>math.exp(i)}
    val total = exponentiated.reduceLeft{_+_} + math.exp(0.0D)
    val ans = new scala.collection.immutable.VectorBuilder[Double]
    ans += math.exp(0.0D)/total 
    ans ++= (0 until d.length).map{i=> exponentiated(i)/total}
    ans.result
  }
  def getOpt(d:IndexedSeq[Double])={
    val logOne=math.log(d.head)
    d.tail.map{i=>math.log(i)/logOne}
  }
  def lower(i:Int)= -10.0
  def upper(i:Int)=10.0
}
case object S extends ParamName{
  def lower(i:Int)=0.0
  def upper(i:Int)=50.0
  def getOpt(d:Seq[Seq[Double]]):IndexedSeq[Double]={
    var initial = d(1)(0) match {
      case 1.0=>d
      case norm => d.map{_.map{_/norm}}
    }
    val ans = for (i<-1 until d.length) yield {
     initial = initial.tail 
     initial.head.take(i)
    }
    ans.flatten
  }
  def getReal(d:Seq[Double])={
    val len = d.length+1
    val matsize = 1 + (math.sqrt(8*len+1).toInt -1)/2 //from triangular root
    val initial = List(1.0).iterator ++ d.iterator
    for (i<-0 until matsize) yield {
      {initial.take(i).toList ++ List.fill(matsize-i)(0.0)}.toIndexedSeq
    }
  }
}

case object BranchLengths extends ParamName{
  def lower(i:Int)=0.0
  def upper(i:Int)=100.0
}

case object Gamma extends ParamName{
  lazy val gamma = new GammaDist
  def lower(i:Int)=0.01
  def upper(i:Int)=1000.0
  def getDist(d:Double,nC:Int)={
    gamma(d,nC)
  }
}


class OptModel[A <: Model](var calc:LikelihoodCalc[A],var tree:Tree,aln:Alignment){
  def m = calc.model
  val myParams:List[(ParamName,Int)] = m.params 

  def update(p:ParamName,d:Double){update(p,d,None)}
  def update(p:ParamName,d:IndexedSeq[Double]){update(p,d,None)}
  def update(p:ParamName,d:IndexedSeq[IndexedSeq[Double]])(implicit m:Manifest[IndexedSeq[IndexedSeq[Double]]]){update(p,d,None)}
  def update(p:ParamName,d:Double,paramIndex:Option[Int]){
    p match {
      case BranchLengths => 
        paramIndex match {
          case None => cantHandle(p,d,paramIndex)
          case Some(i) => tree = tree.setBranchLength(i,d);updatedTree()
        }
      case any => calc=calc.updatedVec(any,Vector(d),paramIndex)
    }
  }

  def update(p:ParamName,value:IndexedSeq[Double],paramIndex:Option[Int]){
    p match {
      case BranchLengths => 
        paramIndex match {
          case None => tree = tree.setBranchLengths(value)
          case Some(i) => cantHandle(p,value,paramIndex)
        }
      case any => calc = calc.updatedVec(p,value,paramIndex)
    }
  }

  def update(p:ParamName,value:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int]=None)(implicit m:Manifest[IndexedSeq[IndexedSeq[Double]]]){
    p match {
      case BranchLengths=> cantHandle(p,value,paramIndex)
      case any => calc = calc.updatedMat(p,value,paramIndex)
    }
  }

  def setOptParam(p:ParamName,value:IndexedSeq[Double],paramIndex:Option[Int]){
    calc = calc.setOptParam(p,value,paramIndex)
  }

  def getOptParam(p:ParamName,paramIndex:Option[Int]=None)={
    m.getOptParam(p,paramIndex)
  }
  def logLikelihood=calc.logLikelihood
  

  def updatedTree(){
    calc = calc updated tree
  }
  def cantHandle(p:ParamName,a:Any,paramIndex:Option[Int]){
    println("Can't handle combination " + p + " " + paramIndex + " " + a)
  }

  def update(t:(ParamName,Int), value:Double){ update(t._1,value,Some(t._2))}
  def update(t:(ParamName,Int), value:IndexedSeq[Double]){ update(t._1,value,Some(t._2))}
  def update(t:(ParamName,Int), value:IndexedSeq[IndexedSeq[Double]])(implicit m:Manifest[IndexedSeq[IndexedSeq[Double]]]){ update(t._1,value,Some(t._2))}

  def apply(t:(ParamName,Int)):Option[IndexedSeq[Double]] = m.getOptParam(t._1,Some(t._2))
  def apply(p:ParamName):Option[IndexedSeq[Double]]=m.getOptParam(p,None)

  def optimiseAll(list:ParamName*){optimiseSeq(list.map{(_,None)})}
  def optimise(list:(ParamName,Option[Int])*){optimiseSeq(list)}
  def optimiseSeq(list:Seq[(ParamName,Option[Int])]){
    val optParams = myParams.filter{t=>
      list.filter{t2=> t2._2.isEmpty || Some(t._2)==t2._2}.map{_._1}.contains(t._1) 
    }
    val s1 = optParams.map{t=>getOptParam(t._1,Some(t._2))}
    if (s1 contains None){error("Not all specified params exist!")}
    val start = s1.map{_.get}.flatten.toArray
    val lengths = s1.map{_.get.length}
    val numArguments = start.length

    import dr.math.{UnivariateFunction,UnivariateMinimum,MultivariateFunction,ConjugateDirectionSearch}
    if (numArguments > 1){
      val func = new MultivariateFunction{
        val lowerBound = optParams.zip(lengths).map{t=> val ((pName,pIndex),len) = t; (0 until len).map{i=> pName.lower(i)}}.flatten
        val upperBound = optParams.zip(lengths).map{t=> val ((pName,pIndex),len) = t; (0 until len).map{i=> pName.upper(i)}}.flatten
        def getLowerBound(i:Int)=lowerBound(i)
        def getUpperBound(i:Int)=upperBound(i)
        val getNumArguments = numArguments
        def evaluate(params:Array[Double])={
          val p2 = params.iterator
          val splitP = lengths.map{l=> p2.take(l).toIndexedSeq}
          optParams.zip(splitP).foreach{ t=> val ((pName,pIndex),values)=t
            setOptParam(pName,values,Some(pIndex))
          }
          val ans = logLikelihood
          println(params.toList + " " + ans)
          -ans
        }
      }
      val search = new ConjugateDirectionSearch
      search.optimize(func,start,1E-4,1E-3)
    }else {
      val func = new UnivariateFunction{
        val param = optParams.head
        val getLowerBound = param._1.lower(0)
        val getUpperBound = param._1.upper(0)
        def evaluate(p:Double)={
          setOptParam(param._1,Vector(p),Some(param._2))
          val ans = logLikelihood
          println(p + " " + ans)
          -ans
        }
      }
      val search = new UnivariateMinimum
      val finalP = search.optimize(start(0),func,1E-4)
      func evaluate finalP
    }
  }
}





