package modiphy.tree
import scala.collection.immutable.Map.WithDefault


abstract class Node
case class Leaf(name:String) extends Node
case class Edge(left:Node,right:Node,dist:Double){
  def hasNode(n:Node) = {left==n || right==n}
  def traverse(n:Node)=n match {
    case `left` => right
    case `right` => left
    case _ => left
  }
  def flip=copy(left=right,right=left)
  def from(n:Node)=n match {
     case `left` => this
     case `right` => flip
  }
  def to(n:Node)=from(n).flip
  def same(that:Edge) = this==that || flip==that
}
class INode extends Node

object TreeTest{
  def main(args:Array[String]){
    val iNode1 = new INode
    val iNode2 = new INode
    val tree = Tree(Edge(Leaf("human"),iNode1,0.10) :: Edge(Leaf("chimp"),iNode1,0.12) :: Edge(iNode1,iNode2,0.21) :: Edge(iNode2,Leaf("gorilla"),0.31) :: Nil)

    println(tree)
    println(tree.treeLength)
  }
}

object Tree{
  def apply(edges:List[Edge],nodes:Set[Node],edgeMap:Map[Node,List[Edge]]):Tree=new Tree(edges,nodes,edgeMap)
  def apply(edges:List[Edge]):Tree={
    val edgeMap = edges.foldLeft(new WithDefault(Map[Node,List[Edge]](),{n:Node=>List[Edge]()})){(m,e)=> m updated (e.left,e::m(e.left)) updated (e.right,e::m(e.right))}
    Tree(edges,edges.foldLeft(Set[Node]()){(s,e)=>(s+e.left)+e.right},edgeMap)
  }
}
class Tree(edges:List[Edge],nodes:Set[Node],edgeMap:Map[Node,List[Edge]]){
  def treeLength = edges.map{_.dist}.foldLeft(0.0D){_+_}
  lazy val iNodes = nodes.filter{_.isInstanceOf[INode]}
  def defaultRoot = iNodes.head
  def getEdges=edgeMap
  def getEdgesTo(n:Node)=edgeMap(n).map{e=>e to n}
  def children(n:Node,dir:Option[Edge]=None)=dir match {
    case None => getEdges(n).map{e=>e from n}
    case Some(ex) => getEdges(n).filter{e=> !(e same ex)}.map{e=>e from n}
  }

  def ancestralTo(n:Node,dir:Edge):Set[Edge]={
    val edgesTo = getEdgesTo(n).filter{e=> !(e same dir)}
    edgesTo.map{e=>ancestralTo(e.left,e)}.foldLeft(edgesTo.toSet){_++_}
  }
    
  def ancestralTo(e:Edge):Set[Edge]=ancestralTo(e.left,e)

  override def toString = toString(defaultRoot)

  def toString(node:Node,direction:Option[Edge]=None):String = { println(node + getEdges(node).map{e=> e from node}.mkString(" "))
   (node,direction) match { 
    case (n:INode,None) => "("+getEdges(node).map{e=>toString(e from n)}.mkString(",")+")"
    case (n:INode,Some(dir)) => "("+getEdges(node).filter{e=> !(e same dir)}.map{e=>toString(e from n)}.mkString(",")+")"
    case (leaf:Leaf,_) => leaf.name
   }}

  def toString(edge:Edge):String = toString(edge.right,Some(edge)) + ":" + edge.dist
}


class Model
abstract class LikelihoodCalculator(t:Tree,m:Model){
  var cache=Map[Edge,Likelihoods]()
  type PartialLikelihoods = Seq[IndexedSeq[Double]]
  type Likelihoods = IndexedSeq[Double]
  type LogLikelihoods = IndexedSeq[Double]
  type Matrix = IndexedSeq[IndexedSeq[Double]]

  def partialLikelihoods(n:Node):PartialLikelihoods={
    combinePartialLikelihoods( (t children n).map{partialLikelihoods} )
  }
  def partialLikelihoods(e:Edge):PartialLikelihoods={
    partialLikelihoodCalc(partialLikelihoods(e.right),matrix(e))
  }
  def logLikelihoods(root:Node,pi:IndexedSeq[Double]):LogLikelihoods={ finalLikelihoods(partialLikelihoods(root),pi).map{math.log} }


  def updated(edge:Edge){cache = cache - edge;(t ancestralTo edge).foreach{e=> cache = cache - e}}
  def updated{cache=Map[Edge,Likelihoods]()}

  def matrix(edge:Edge):Matrix
  def partialLikelihoodCalc(end:PartialLikelihoods,matrix:Matrix):PartialLikelihoods    
  def combinePartialLikelihoods(intermediates:List[PartialLikelihoods]):PartialLikelihoods
  def finalLikelihoods(partial:PartialLikelihoods,pi:IndexedSeq[Double]):Likelihoods
}

trait ColtLikelihoodCalc{
  import cern.colt.matrix._
  val func = new cern.colt.function.DoubleDoubleFunction{def apply(x:Double,y:Double)=x*y}
  val fact1D = cern.colt.matrix.DoubleFactory1D.dense
  val fact2D = cern.colt.matrix.DoubleFactory2D.dense


  def partialLikelihoodCalc(end:List[DoubleMatrix1D],matrix:DoubleMatrix2D)={
    val width = matrix.rows
    val rows = new Array[DoubleMatrix1D](width)
    (0 until width).foreach{i=> rows(i)=matrix viewRow i}
    end.map{siteVector=>
        val ret = fact1D.make(width)

        for (i<-0 until width){
          ret set (i,siteVector.zDotProduct(rows(i)))
        }
        ret
      }
  }
  def combinePartialLikelihoods(intermediates:List[List[DoubleMatrix1D]])={
    val ans = intermediates.head
    intermediates.tail.foreach{list2=>
      ans.zip(list2).map{t=> // not really a map but used for parallel reasons
        val (vec,vec2)=t
        vec.assign(vec2,func)
      }
    }
    ans

  }

  def likelihoods(partial:List[DoubleMatrix1D],pi:Seq[Double])={
    partial.map{vec=>
      val ans = pi.iterator.zipWithIndex.map{t=>
        val(p,i)=t
        vec.get(i)*p
      }.foldLeft(0.0D){_+_}
      ans
    }
  }
}
