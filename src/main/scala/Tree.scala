package modiphy.tree
import scala.collection.immutable.Map.WithDefault
import modiphy.alignment._
import modiphy.model._
import modiphy.alignment.GlobalAlphabet._
import scala.collection.LinearSeq

import modiphy.opt._

object Parallel{
  var on = true
  import jsr166y._
  lazy val forkJoinPool = new ForkJoinPool
  lazy val threshold = -1
}
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
  def different(that:Edge)= !(same(that))
}
class INode extends Node

object TreeTest{
  def main(args:Array[String]){
    import modiphy.test.ModelData._
    import modiphy.math.constants._
    /*
    val iNode1 = new INode
    val iNode2 = new INode
    val tree = Tree(Edge(Leaf("human"),iNode1,0.10) :: Edge(Leaf("chimp"),iNode1,0.12) :: Edge(iNode1,iNode2,0.21) :: Edge(iNode2,Leaf("gorilla"),0.31) :: Nil)

    println(tree)
    println(tree.treeLength)
    */
    /*

    if(args.length >1 && (args(1) startsWith "par")){Parallel.on=true}
    if(args.length >1 && (args(1) startsWith "nopar")){Parallel.on=false}
    val tree = Tree(treeStr)
    val aln = Fasta(alnStr).parseWith(AminoAcid)
    */
/*
    val model = BasicLikelihoodModel(WAG.pi,WAG.S)


    (0 until args(0).toInt).foreach{i=>
      val lkl = new SimpleLikelihoodCalc(tree,model,aln) 
      println(lkl.logLikelihood)
    }
    */
    val tree = Tree(tufaTree)
    val aln  = new Fasta(tufaAln.lines) parseWith AminoAcid
 
    val model = GammaModel(aln.frequencies,WAG.S,0.5,4)
    val lkl = new MixtureLikelihoodCalc(tree,aln,model)
    val optModel = new OptModel(lkl,tree,aln)
    optModel optimiseAll Gamma
    optModel optimiseAll (Pi,Gamma)
    println(optModel)

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
  import java.io.File
  def apply(file:File):Tree = {
    apply(scala.io.Source.fromFile(file).getLines.map{_.trim}.mkString(""))
  }
}

trait TreePosition{
  val get:Node
  def neighbours:Seq[TreePosition]
}
trait RootedTreePosition extends TreePosition{
  import modiphy.util.Memo
  val children:Seq[TreePositionDir]
  def neighbours:Seq[RootedTreePosition] = children
  val ancestralTo:Memo[Node,Boolean] = Memo[Node,Boolean]({n=>
    neighbours.foldLeft(false){(bool,neighbour)=>
      bool || n==neighbour.get || neighbour.ancestralTo(n) //this last one could blow the stack....
    } 
  })
  def hashCodeFunc:Int = children.map{_.hashCode}.foldLeft(get.hashCode + 41){(h1,child)=>
    h1 * 41 + child.hashCode
  }
  override lazy val hashCode = hashCodeFunc
  override def equals(o:Any)={
    o match {
      case that:TreePositionDir=>false
      case that:RootedTreePosition=>that.get==get && (that.tree eq tree)
      case _ => false
    }
  }
  val tree:Tree
  val leaves:Seq[Leaf]
  val id:Int
}
trait TreePositionDir extends RootedTreePosition{
  val upEdge:Edge
  override lazy val hashCode:Int = super.hashCodeFunc + 41 * upEdge.hashCode
  override def equals(o:Any)={
    o match {
      case that:TreePositionDir=>that.upEdge==upEdge && that.get==get && (that.tree eq tree)
      case _ => false
    }
  }
  val tree:Tree
}

class Tree(val edges:IndexedSeq[Edge],
  nodes:Set[Node],
  edgeMap:Map[Node,List[Edge]],
  startiNodes:Option[Set[INode]] = None,
  root:Option[INode]=None,
  startLabels:Option[Map[String,Node]] = None){

  lazy val branchLength = edges.map{_.dist}.toIndexedSeq


  def drop(leafName:String)={
    val hasEdge = edgeMap.filter{t=> t._1==Leaf(leafName)}//.head._2.head
    if (hasEdge.isEmpty){
      this
    }else {
      val upEdge = hasEdge.head._2.head
      val iNode = upEdge from Leaf(leafName) right

      val newNodes = nodes.filter{n=> n!=Leaf(leafName) } 
      val newEdges = edgeMap.filter{t=> t._1!=Leaf(leafName)}.map{t=> (t._1, t._2.filter{e=>e different upEdge})}
      new Tree(edges.filter{e => e different upEdge},newNodes,newEdges,None,root,None) dropIfSafe (iNode.asInstanceOf[INode])
    }
  }
  def dropIfSafe(iNode:INode)={
    if (edgeMap(iNode).length >2 ){
      this
    }else {
      val uselessEdges = edgeMap(iNode)
      val len = uselessEdges.map{_.dist}.reduceLeft{_+_}
      val ends = uselessEdges.map{_ from iNode right}

      val newEdge = Edge(ends(0),ends(1),len)
        val newEdgeMap = edgeMap.filter{t=> t._1 != iNode}.updated(ends(0), newEdge::edgeMap(ends(0)).filter{e=> e different uselessEdges(0)}).updated(ends(1),newEdge::edgeMap(ends(1)).filter{e=> e different uselessEdges(1)})

        new Tree(edges.filter{e=>uselessEdges.find(e2=> e same e2).isEmpty} :+ newEdge, nodes.filter{n=> n!=iNode},newEdgeMap) 
    }
  }

  def restrictTo(list:Seq[String])={
    val remove = leafNodes.map{_.name}.filter{name => !list.contains{name}}
    remove.foldLeft(this){(tree,leaf)=> tree.drop(leaf)}
  }
  

  def numNodes = nodes.size
  def numLeaves = nodes.filter{n => n.isInstanceOf[Leaf]}.size
  

  val iNodes = startiNodes.getOrElse(nodes.filter{_.isInstanceOf[INode]}.map{_.asInstanceOf[INode]}).toIndexedSeq
  def nodeID(n:INode)=iNodes.indexOf(n)
  def edgeID(e:Edge)={
    edges.indexOf(e) match {
      case -1 => edges.indexOf(e.flip) + edges.length
      case i => i
    }
  }


  val leafNodes = nodes.filter{_.isInstanceOf[Leaf]}.map{_.asInstanceOf[Leaf]}
  val defaultRoot = root getOrElse iNodes.head
  val labels = startLabels getOrElse nodes.foldLeft(Map[String,Node]()){(m,n)=> if (n.id.isDefined){m updated (n.id.get,n)}else{m}}
  lazy val edgeSet = edges.foldLeft(Set[Edge]()){_+_}

  def hasEdge(e:Edge)=edgeSet.contains(e)

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
    val parent =this
    new TreePositionDir{
      val get = n
      lazy val children = edgeMap(n).filter{e=> ! (e same dir)}.map{e=>traverseDown(e from n right,e)}
      val upEdge = dir
      val leaves = parent.leaves(n,dir)
      val tree = parent
      val id = edgeID(dir to n)
    }
  }

  def traverseDown(n:INode):RootedTreePosition={
    val parent =this
    new RootedTreePosition{
      val get = n
      lazy val children = edgeMap(n).map{e=>traverseDown(e from n right,e)}
      val leaves = parent.leaves(n)
      val tree = parent
      val id = nodeID(n) + edges.length * 2
    }
  }

  val maxID = edges.length * 2 + iNodes.length

  def copy(
    newEdges:IndexedSeq[Edge] = edges,
    newNode:Set[Node] = nodes,
    newEdgeMap:Map[Node,List[Edge]] = edgeMap,
    newINodes:Option[Set[INode]] = Some(iNodes.toSet),
    root:Option[INode] = Some(defaultRoot)
  ) = new Tree(newEdges,newNode,newEdgeMap,newINodes,root)
  def reRoot(n:INode)=copy(root=Some(n))
  def setBranchLength(i:Int,d:Double)={
    val oldEdge = edges(i)
    val e = oldEdge
    val newEdge = edges(i).copy(dist=d)
    val myNewEdgeMap = edgeMap.updated(e.left,newEdge::edgeMap(e.left).filter{_!=oldEdge}).updated(e.right,newEdge::edgeMap(e.right).filter{_!=oldEdge})
    val nE = edges.updated(i,newEdge)
    assert(nE.length==edges.length)
      val ans = copy(newEdges = edges.updated(i,newEdge),newEdgeMap = myNewEdgeMap)
  //    println("New tree " + i + " " + d + " " + ans)
    ans
  }
  def setBranchLengths(vec:Seq[Double])={
    var myNewEdgeMap = edgeMap
    if (vec.length < edges.length){error("Tree with " + edges.length + " edges but only " + vec.length + " lengths specified")}
    val edges2 = edges.zip(vec).map{ t=>val (e,d)=t
      if (e.dist!=d){
        val ans = e.copy(dist=d)
        myNewEdgeMap = myNewEdgeMap.updated(e.left,ans::myNewEdgeMap(e.left).filter{_!=e}).updated(e.right,ans::myNewEdgeMap(e.right).filter{_!=e})
        ans
      }else {
        e
      }
    }
    assert(edges2.length==edges.length)
    copy(newEdges=edges2,newEdgeMap = myNewEdgeMap)
  }
  lazy val getBranchLengths = edges.map{_.dist}

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

import scalaz.Scalaz._
  lazy val leafCache:((Node,Option[Edge]))=>IndexedSeq[Leaf] = mutableHashMapMemo{t:(Node,Option[Edge])=>
    t match {
      case (l:Leaf,_)=>Vector(l)
      case (n:Node,e)=>children(n,e).map{e=> leafCache((e from n right,Some(e)))}.foldLeft(Vector[Leaf]()){_++_}
    }
  }
  def leaves(n:Node,e:Edge):Seq[Leaf]=this.synchronized{leafCache(n,Some(e))}
  def leaves(n:Node):Seq[Leaf]=this.synchronized{leafCache(n,None)}
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

  def differences(t2:Tree):Seq[Edge]={
    if (t2 eq this){
      List()
    }else {
      edges.filter{e=> ! (t2 hasEdge e)}
    }
  }

}

