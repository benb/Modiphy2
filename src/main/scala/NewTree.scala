package modiphy.tree
import scala.collection.immutable.Map.WithDefault
import modiphy.alignment._
import modiphy.model._
import modiphy.alignment.GlobalAlphabet._
import scala.collection.LinearSeq

import modiphy.opt._

object Parallel{
  val on = false
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
}
class INode extends Node

object TreeTest{
  def main(args:Array[String]){
    /*
    val iNode1 = new INode
    val iNode2 = new INode
    val tree = Tree(Edge(Leaf("human"),iNode1,0.10) :: Edge(Leaf("chimp"),iNode1,0.12) :: Edge(iNode1,iNode2,0.21) :: Edge(iNode2,Leaf("gorilla"),0.31) :: Nil)

    println(tree)
    println(tree.treeLength)
    */
    import modiphy.test.ModelData._
    import modiphy.math.constants._

    val tree = Tree(treeStr)
    val aln = Fasta(alnStr).parseWith(AminoAcid)
    val model = BasicLikelihoodModel(WAG.pi,WAG.S)

    (0 until args(0).toInt).foreach{i=>
      val lkl = new SimpleLikelihoodCalc(tree,model,aln) 
      println(lkl.logLikelihood)
    }
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

class Tree(edges:IndexedSeq[Edge],
  nodes:Set[Node],
  edgeMap:Map[Node,List[Edge]],
  startiNodes:Option[Set[INode]] = None,
  root:Option[INode]=None,
  startLabels:Option[Map[String,Node]] = None){

  lazy val branchLength = edges.map{_.dist}.toIndexedSeq

  

  

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
      copy(newEdges = edges.updated(i,newEdge),newEdgeMap = myNewEdgeMap)
  }
  def setBranchLengths(vec:Seq[Double])={
    var myNewEdgeMap = edgeMap
    val edges2 = edges.zip(vec).map{ t=>val (e,d)=t
      if (e.dist!=d){
        val ans = e.copy(dist=d)
        myNewEdgeMap = myNewEdgeMap.updated(e.left,ans::myNewEdgeMap(e.left).filter{_!=e}).updated(e.right,ans::edgeMap(e.right).filter{_!=e})
        ans
      }else {
        e
      }
    }
    copy(newEdges=edges2,newEdgeMap = myNewEdgeMap)
  }
  lazy val getBranchLength = edges.map{_.dist}

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

import modiphy.util.Memo
  lazy val leafCache:Memo[(Node,Option[Edge]),Seq[Leaf]]=Memo[(Node,Option[Edge]),Seq[Leaf]]({t=>
    t match {
      case (l:Leaf,_)=>Vector(l)
      case (n:Node,e)=>children(n,e).map{e=> leafCache((e from n right,Some(e)))}.foldLeft(Vector[Leaf]()){_++_}
    }
  })
  def leaves(n:Node,e:Edge):Seq[Leaf]=leafCache(n,Some(e))
  def leaves(n:Node):Seq[Leaf]=leafCache(n,None)
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



object LikelihoodTypes{
  type Pattern=String=>Letter
  type PartialLikelihoods = LinearSeq[Double]
  type Likelihood = Double
  type LogLikelihood = Double
  type Matrix = LinearSeq[LinearSeq[Double]]
}
import LikelihoodTypes._

object SimpleLikelihoodCalc{
  type Cache=Map[RootedTreePosition,Map[Seq[Letter],PartialLikelihoods]]
}

trait LikelihoodFactory{
  def apply:LikelihoodEngine
}
object DefaultLikelihoodFactory{
  var default:LikelihoodFactory = IndexedSeqLikelihoodFactory
  def apply = default.apply
  def setDefault(lkl:LikelihoodFactory){default=lkl}
}
object ColtLikelihoodFactory extends LikelihoodFactory{
  def apply=new ColtLikelihoodCalc
}
object IndexedSeqLikelihoodFactory extends LikelihoodFactory{
  def apply=new IndexedSeqLikelihoodCalc
}


class MixtureLikelihoodCalc(priors:Seq[Double],tree:Tree,aln:Alignment,m:StdMixtureModel,lkl:Option[Seq[SimpleLikelihoodCalc]]=None) extends LikelihoodCalc[StdMixtureModel]{
  import scala.actors.Futures.future
  val models = m.models
  val lklCalc = lkl.getOrElse{models.map{new SimpleLikelihoodCalc(tree,_,aln)}}
  
  lazy val logLikelihood={
    val likelihoods = lklCalc.zip(priors).map{t=> t._1.likelihoods().map{_ * t._2}}.map{_.iterator}
    var ans =  0.0
    val patternCount = aln.countList.iterator
    while (likelihoods.head.hasNext){
      val row = likelihoods.map{_.next}
      ans = ans + math.log(row.reduceLeft{_+_}) * patternCount.next
    }
    ans
  }
  def updated(t:Tree)=new MixtureLikelihoodCalc(priors,t,aln,m,lkl)
  def updated(m:StdMixtureModel)=new MixtureLikelihoodCalc(priors,tree,aln,m,lkl)

   def updatedVec(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int])={ updated(m.updatedVec(p,vec,paramIndex)) }
   def updatedMat(p:ParamName,mat:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int])={ updated(m.updatedMat(p,mat,paramIndex)) }
}

/*
trait CachingLikelihoodCalc[A <: Model] extends LikelihoodCalc[A]{
   val partialLikelihoodCache:modiphy.util.Memo[RootedTreePosition,LinearSeq[PartialLikelihoods]]=new modiphy.util.ArrayMemo[RootedTreePosition,LinearSeq[PartialLikelihoods]]({treePos=>treePos.id},tree.maxID)({ treePos => realPartialLikelihoodCalc(treePos)})
   override def partialLikelihoods(treePos:RootedTreePosition)=partialLikelihoodCache(treePos)  
   def tree:Tree
   def realPartialLikelihoodCalc(treePos:RootedTreePosition):LinearSeq[PartialLikelihoods]
}*/
trait LikelihoodCalc[A <: Model]{
//   def partialLikelihoods(treePos:RootedTreePosition):LinearSeq[PartialLikelihoods]
   def updated(t:Tree):LikelihoodCalc[A]
   def updated(m:A):LikelihoodCalc[A]
   def logLikelihood:Double
   def updatedVec(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int]):LikelihoodCalc[A]
   def updatedMat(p:ParamName,mat:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int]):LikelihoodCalc[A]
}
class SimpleLikelihoodCalc(val tree:Tree,m:StdModel,val aln:Alignment,val engine:LikelihoodEngine=DefaultLikelihoodFactory.apply) extends LikelihoodCalc[StdModel]{
  import SimpleLikelihoodCalc._
  
  import engine.combinePartialLikelihoods
  import engine.partialLikelihoodCalc
  import engine.finalLikelihood

   def updatedVec(p:ParamName,vec:IndexedSeq[Double],paramIndex:Option[Int])={ updated(m.updatedVec(p,vec,paramIndex)) }
   def updatedMat(p:ParamName,mat:IndexedSeq[IndexedSeq[Double]],paramIndex:Option[Int])={ updated(m.updatedMat(p,mat,paramIndex)) }
  def realPartialLikelihoodCalc(treePos:RootedTreePosition):LinearSeq[PartialLikelihoods]={
    val p = aln.patternList
       treePos.get match {
          case n:INode=>
            combinePartialLikelihoods(
              treePos.children.toList.map{tp=>
              val plStart:LinearSeq[PartialLikelihoods] = partialLikelihoods(tp)
              partialLikelihoodCalc(plStart,m(tp.upEdge)) 
            }
          )
          case l:Leaf=>p.map{p2=>leafPartialLikelihoods(p2(l.id.get))}
        }
    }
    import scalaz._
    import Scalaz._

//   val partialLikelihoods = mutableHashMapMemo(realPartialLikelihoodCalc)
//   val partialLikelihoods:modiphy.util.Memo[RootedTreePosition,LinearSeq[PartialLikelihoods]]=modiphy.util.Memo[RootedTreePosition,LinearSeq[PartialLikelihoods]]({ treePos => realPartialLikelihoodCalc(treePos)})
//   val partialLikelihoods:modiphy.util.Memo[RootedTreePosition,LinearSeq[PartialLikelihoods]]=new modiphy.util.ArrayMemo[RootedTreePosition,LinearSeq[PartialLikelihoods]]({treePos=>treePos.id},tree.maxID)({ treePos => realPartialLikelihoodCalc(treePos)})
    def partialLikelihoods(treePos:RootedTreePosition)=realPartialLikelihoodCalc(treePos)  


  def likelihoods(root:RootedTreePosition=tree.traverseDown(tree.defaultRoot)):Seq[Double]={
    finalLikelihood(partialLikelihoods(root),m.pi(root.get))
  }

  def logLikelihoodRoot(root:RootedTreePosition=tree.traverseDown(tree.defaultRoot)):Double={
    likelihoods(root).zip(aln.countList).map{t=>math.log(t._1)*t._2}.reduceLeft{_+_}
  }
  lazy val logLikelihood:Double=logLikelihoodRoot()

  val leafPartialLikelihoods=immutableHashMapMemo{l:Letter=>l match {
    case a if (a.isReal) => List.fill(l.alphabet.length)(0.0).updated(l.id,1.0)
    case a => List.fill(l.alphabet.length)(1.0)
  }}


  def factory(t:modiphy.tree.Tree,m:StdModel,aln:Alignment) = new SimpleLikelihoodCalc(t,m,aln)

  def updated(t:modiphy.tree.Tree)={
    /*
    val updatedCache = tree.differences(t).foldLeft(cache){(c2,e)=>
        c2.keys.foldLeft(c2){(map,tp)=> 
         if (tp ancestralTo e.left){ // if I am ancestral to edge e then both nodes are my descendents
           map - tp
         }else {
           map
         }
      }
    }*/
    factory(t,m,aln)
  }

  def updated(m:StdModel)={
    factory(tree,m,aln)//,Map[RootedTreePosition,Map[Seq[Letter],PartialLikelihoods]]())
  }
}

trait LikelihoodEngine{
  def combinePartialLikelihoods(intermediates:LinearSeq[LinearSeq[PartialLikelihoods]]):LinearSeq[PartialLikelihoods]
  def partialLikelihoodCalc(end:LinearSeq[PartialLikelihoods],matrix:Matrix):LinearSeq[PartialLikelihoods]    
  def finalLikelihood(partial:LinearSeq[PartialLikelihoods],pi:IndexedSeq[Double]):LinearSeq[Likelihood]
}
class ColtLikelihoodCalc extends LikelihoodEngine{
  import SimpleLikelihoodCalc.Cache
  import cern.colt.matrix._
  val func = new cern.colt.function.DoubleDoubleFunction{def apply(x:Double,y:Double)=x*y}

  import modiphy.math.EnhancedMatrix._

  def combinePartialLikelihoods(intermediates:LinearSeq[LinearSeq[PartialLikelihoods]]):LinearSeq[PartialLikelihoods] = {
    val myInter:LinearSeq[LinearSeq[DoubleMatrix1D]] = intermediates.map{list => 
      list.map{p:PartialLikelihoods=>
        Seq2Vec(p)
      }
    }
    coltCombinePartialLikelihoods(myInter).map{Vec2Seq}
  }
  def partialLikelihoodCalc(end:LinearSeq[PartialLikelihoods],matrix:Matrix)={
    coltPartialLikelihoodCalc(end.map(Seq2Vec),matrix).map{Vec2Seq}
  }
  def finalLikelihood(partial:LinearSeq[PartialLikelihoods],pi:IndexedSeq[Double])={
    coltLikelihoods(partial.map(Seq2Vec),pi)
  }


  def coltPartialLikelihoodCalc(end:LinearSeq[DoubleMatrix1D],matrix:DoubleMatrix2D)={
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
  def coltCombinePartialLikelihoods(intermediates:LinearSeq[LinearSeq[DoubleMatrix1D]])={
    val ans = intermediates.head
    intermediates.tail.foreach{list2=>
      ans.zip(list2).map{t=> // not really a map but used for parallel reasons
        val (vec,vec2)=t
        vec.assign(vec2,func)
      }
    }
    ans
  }

  def coltLikelihoods(partial:LinearSeq[DoubleMatrix1D],pi:Seq[Double])={
    partial.map{vec=>
      val ans = pi.iterator.zipWithIndex.map{t=>
        val(p,i)=t
        vec.get(i)*p
      }.foldLeft(0.0D){_+_}
      ans
    }
  }
}

class IndexedSeqLikelihoodCalc extends LikelihoodEngine{
  import SimpleLikelihoodCalc.Cache

  import modiphy.math.EnhancedMatrix._
  import modiphy.math.EnhancedMatrix

  def combinePartialLikelihoods(intermediates:LinearSeq[LinearSeq[PartialLikelihoods]]):LinearSeq[PartialLikelihoods] = {
    //TODO opt
    var ans = intermediates.head
    intermediates.tail.foreach{list2=>
      ans = ans.zip(list2).map{t=>
        t._1 prod t._2
      }
    }
    ans
  }
  
  def partialLikelihoodCalc(end:LinearSeq[PartialLikelihoods],matrix:Matrix)={
    //TODO opt
    new EnhancedMatrix(end) multipleDotProduct matrix
  }

  def finalLikelihood(partial:LinearSeq[PartialLikelihoods],pi:IndexedSeq[Double])={
    val piList = pi.toList
    partial.map{vec=>
      piList.dotProduct(vec)
    }
  }
}
