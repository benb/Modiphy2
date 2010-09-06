package modiphy.opt
import modiphy.tree._
import modiphy.alignment._
import modiphy.model._

trait ParamName{
  def apply(i:Int)=(this,i)
  def lower(i:Int):Double
  def upper(i:Int):Double 
}

case object Pi extends ParamName{
  def getReal(d:IndexedSeq[Double])={
    val exponentiated =  d.map{i=>Math.exp(i)}
    val total = exponentiated.reduceLeft{_+_} + Math.exp(0.0D)
    val ans = new scala.collection.immutable.VectorBuilder[Double]
    ans += Math.exp(0.0D)/total 
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
  def upper(i:Int)=1000.0
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
  def upper(i:Int)=50.0
}
abstract class OptModel(var m:StdModel,var calc:LikelihoodCalc,var tree:Tree,aln:Alignment){
  val myParams:List[(ParamName,Int)] = m.params 
  def update(p:ParamName,value:Any,paramIndex:Option[Int]=None){
    p match {
      case BranchLengths => 
        paramIndex match {
          case None => value match {
            case v:IndexedSeq[Double] => tree = tree.setBranchLengths(v);updatedTree()
            case a => cantHandle(p,a,paramIndex)
          }
          case Some(i) => value match {
            case d:Double => tree = tree.setBranchLength(i,d);updatedTree()
            case a => cantHandle(p,a,paramIndex)
          }
        }
      case any => value match {
        case d:Double => m = m.updatedVec(p,Vector(d),paramIndex);updatedModel
        case vec:IndexedSeq[Double] => m = m updatedVec(p,vec,paramIndex);updatedModel
        case mat:IndexedSeq[IndexedSeq[Double]] => m = m updatedMat(p,mat,paramIndex);updatedModel
        case vec:Seq[Double] => update(p,vec.toIndexedSeq,paramIndex)
        case any => cantHandle(p,value,paramIndex)
      }
    }
  }
  def logLikelihood=calc.logLikelihood
  

  def updatedTree(){
    calc = calc updated tree
  }
  def updatedModel(){
    calc = calc updated m
  }
  def cantHandle(p:ParamName,a:Any,paramIndex:Option[Int]){
    println("Can't handle combination " + p + " " + paramIndex + " " + a)
  }

  def update(t:(ParamName,Int), value:Any){
    update(t._1,value,Some(t._2))
  }
  def apply(t:(ParamName,Int)):IndexedSeq[Double] = m.getOptParam(t._1,Some(t._2))
  def apply(p:ParamName):IndexedSeq[Double]=m.getOptParam(p,None)
}

