package modiphy.tree
import scala.collection.immutable.Map.WithDefault
import modiphy.util.Memo
import modiphy.alignment._
import modiphy.alignment.GlobalAlphabet._


abstract class Node{
 def id:Option[String]=None
}
case class Leaf(name:String) extends Node{override val id=Some(name)}
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
  def apply(edges:IndexedSeq[Edge],nodes:Set[Node],edgeMap:Map[Node,List[Edge]]):Tree=new Tree(edges,nodes,edgeMap)
  def apply(edges:List[Edge]):Tree={
    val edgeMap = edges.foldLeft(new WithDefault(Map[Node,List[Edge]](),{n:Node=>List[Edge]()})){(m,e)=> m updated (e.left,e::m(e.left)) updated (e.right,e::m(e.right))}
    Tree(edges.toIndexedSeq,edges.foldLeft(Set[Node]()){(s,e)=>(s+e.left)+e.right},edgeMap)
  }
  def apply(newick:String):Tree = {
    new TreeParser{def parseAll=parse(tree,newick)}.parseAll.get
  }
}

trait TreePosition{
  val get:Node
  val neighbours:Seq[TreePosition]
}
trait RootedTreePosition extends TreePosition{
  val children:Seq[TreePositionDir]
  val neighbours:Seq[TreePosition] = children
}
trait TreePositionDir extends RootedTreePosition{
  val upEdge:Edge
}

class Tree(edges:IndexedSeq[Edge],
  nodes:Set[Node],
  edgeMap:Map[Node,List[Edge]],
  startiNodes:Option[Set[INode]] = None,
  root:Option[Node]=None,
  startLabels:Option[Map[String,Node]] = None){

  lazy val branchLength = edges.map{_.dist}.toIndexedSeq


  

  val iNodes = startiNodes.getOrElse(nodes.filter{_.isInstanceOf[INode]}.map{_.asInstanceOf[INode]})
  val leafNodes = nodes.filter{_.isInstanceOf[Leaf]}.map{_.asInstanceOf[Leaf]}
  def defaultRoot = root getOrElse iNodes.head
  val labels = startLabels getOrElse nodes.foldLeft(Map[String,Node]()){(m,n)=> if (n.id.isDefined){m updated (n.id.get,n)}else{m}}

  def traverseFrom(n:Node):Option[TreePosition]=Some(new TreePosition{
    val get = n
    lazy val neighbours = edgeMap(n).map{e=>traverseFrom(e from n right).get}
  })

  def traverseFrom(s:String):Option[TreePosition]={
    if (labels contains s){
      traverseFrom(labels(s))
    }else {
      None
    }
  }

  def traverseDown(n:Node,dir:Edge):TreePositionDir={
    new TreePositionDir{
      val get = n
      lazy val children = edgeMap(n).filter{e=> ! (e same dir)}.map{e=>traverseDown(e from n right,e)}
      val upEdge = dir
    }
  }


  def copy(
    newEdges:IndexedSeq[Edge] = edges,
    newNode:Set[Node] = nodes,
    newEdgeMap:Map[Node,List[Edge]] = edgeMap,
    newINodes:Option[Set[INode]] = Some(iNodes),
    root:Option[Node] = Some(defaultRoot)
  ) = new Tree(newEdges,newNode,newEdgeMap,newINodes,root)
  def reRoot(n:Node)=copy(root=Some(n))
  def setBranchLength(i:Int,d:Double)={
    copy(newEdges = edges.updated(i,edges(i).copy(dist=d)))
  }

  def treeLength = edges.map{_.dist}.foldLeft(0.0D){_+_}
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

  override def toString = toString(defaultRoot) + ";"

  lazy val leaves:Memo[(Node,Edge),Seq[Leaf]]=new Memo[(Node,Edge),Seq[Leaf]]({t=>val (n,e)=t
    n match {
      case l:Leaf=>Vector(l)
      case n:Node=>traverseDown(n,e).children.map{tp:TreePositionDir=> leaves((tp.get,tp.upEdge))}.foldLeft(Vector[Leaf]()){_++_}
    }
  })
  def leaves(treePos:TreePositionDir):Seq[Leaf]=leaves(treePos.get,treePos.upEdge)
  def leaves(treePos:RootedTreePosition):Seq[Leaf]=
    treePos match {
      case tP:TreePositionDir => leaves(tP.get,tP.upEdge)
      case root:RootedTreePosition => leafNodes.toList
    }

  def toString(node:Node,direction:Option[Edge]=None):String = { 
   (node,direction) match { 
    case (n:INode,None) => "("+getEdges(node).map{e=>toString(e from n)}.mkString(",")+")"
    case (n:INode,Some(dir)) => "("+getEdges(node).filter{e=> !(e same dir)}.map{e=>toString(e from n)}.mkString(",")+")"
    case (leaf:Leaf,_) => leaf.name
   }}

  def toString(edge:Edge):String = toString(edge.right,Some(edge)) + ":" + edge.dist

}


trait Model{
  def apply(e:Edge):IndexedSeq[IndexedSeq[Double]]
  def pi(n:Node):IndexedSeq[Double]
}

object LikelihoodTypes{
  type Pattern=Leaf=>Letter
  type PartialLikelihoods = IndexedSeq[Double]
  type Likelihood = Double
  type LogLikelihood = Double
  type Matrix = IndexedSeq[IndexedSeq[Double]]
}
import LikelihoodTypes._

abstract class SimpleLikelihoodCalc(tree:Tree,m:Model){
  
  var cache = Map[(RootedTreePosition,Seq[Letter]),PartialLikelihoods]()


  def partialLikelihoods(treePos:RootedTreePosition,p:Pattern):PartialLikelihoods={
    val myPatterns = tree.leaves(treePos).map{p}
    if (cache contains ((treePos,myPatterns))){
      cache((treePos,myPatterns))
    }else {
    val ans = treePos.get match {
      case n:INode=>
        combinePartialLikelihoods(
          treePos.children.toList.map{tp=>
          val plStart = partialLikelihoods(tp,p)
          //calculate PL along branch
          partialLikelihoodCalc(plStart,m(tp.upEdge)) 
          }
        )
      case l:Leaf=>leafPartialLikelihoods(p(l))
    }
        
    cache = cache updated ((treePos,myPatterns),ans)
    ans
    }
  }

  def likelihood(root:RootedTreePosition,p:Pattern):Likelihood={
     finalLikelihood(partialLikelihoods(root,p),m.pi(root.get))
  }
  def leafPartialLikelihoods(l:Letter):PartialLikelihoods = {
    val empty = Vector.fill(l.alphabet.length)(0.0)
    empty updated (l.id,1.0)
  }

  def combinePartialLikelihoods(intermediates:List[PartialLikelihoods]):PartialLikelihoods
  def partialLikelihoodCalc(end:PartialLikelihoods,matrix:Matrix):PartialLikelihoods    
  def finalLikelihood(partial:PartialLikelihoods,pi:IndexedSeq[Double]):Likelihood
}

trait LikelihoodEngine{
  def combinePartialLikelihoods(intermediates:List[PartialLikelihoods]):PartialLikelihoods
  def partialLikelihoodCalc(end:PartialLikelihoods,matrix:Matrix):PartialLikelihoods    
  def finalLikelihood(partial:PartialLikelihoods,pi:IndexedSeq[Double]):Likelihood
}
trait ColtLikelihoodCalc extends LikelihoodEngine{
  import cern.colt.matrix._
  val func = new cern.colt.function.DoubleDoubleFunction{def apply(x:Double,y:Double)=x*y}
  val fact1D = cern.colt.matrix.DoubleFactory1D.dense
  val fact2D = cern.colt.matrix.DoubleFactory2D.dense

  implicit def Seq2Vec(s:IndexedSeq[Double])=fact1D.make(s.toArray)
  implicit def Vec2Seq(v:DoubleMatrix1D)=v.toArray.toIndexedSeq
  implicit def Seq2Mat(s:IndexedSeq[IndexedSeq[Double]])=fact2D.make(s.map{_.toArray}.toArray)
  implicit def Mat2Seq(m:DoubleMatrix2D)=m.toArray.map{_.toIndexedSeq}.toIndexedSeq

 
  def combinePartialLikelihoods(intermediates:List[PartialLikelihoods]):PartialLikelihoods = {
    val myInter = intermediates.map{p:PartialLikelihoods=>List(Seq2Vec(p))}
    coltCombinePartialLikelihoods(myInter).head
  }
  def partialLikelihoodCalc(end:PartialLikelihoods,matrix:Matrix):PartialLikelihoods={
    coltPartialLikelihoodCalc(List(Seq2Vec(end)),matrix).head
  }
  def finalLikelihood(partial:PartialLikelihoods,pi:IndexedSeq[Double]):Likelihood={
    coltLikelihoods(List(Seq2Vec(partial)),pi).head
  }


  def coltPartialLikelihoodCalc(end:List[DoubleMatrix1D],matrix:DoubleMatrix2D)={
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
  def coltCombinePartialLikelihoods(intermediates:List[List[DoubleMatrix1D]])={
    val ans = intermediates.head
    intermediates.tail.foreach{list2=>
      ans.zip(list2).map{t=> // not really a map but used for parallel reasons
        val (vec,vec2)=t
        vec.assign(vec2,func)
      }
    }
    ans

  }

  def coltLikelihoods(partial:List[DoubleMatrix1D],pi:Seq[Double])={
    partial.map{vec=>
      val ans = pi.iterator.zipWithIndex.map{t=>
        val(p,i)=t
        vec.get(i)*p
      }.foldLeft(0.0D){_+_}
      ans
    }
  }
}
