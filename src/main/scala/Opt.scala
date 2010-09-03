package modiphy.opt
import modiphy.tree._
import modiphy.alignment._
import modiphy.model._

trait ParamName{
  def apply(i:Int)=(this,i)
}

case object Pi extends ParamName
case object S extends ParamName
case object BranchLengths extends ParamName

abstract class OptModel(m:Model,tree:Tree,aln:Alignment){
  val myParams:List[(ParamName,Int)] = m.params 
  def update(p:ParamName,value:Any){
     
  }
  def update(t:(ParamName,Int), value:Any){
  }
  def apply(t:(ParamName,Int))
}

