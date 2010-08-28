package modiphy.tree
import scala.collection.immutable.Map.WithDefault
import modiphy.util.Memo
import modiphy.alignment._
import modiphy.model._
import modiphy.alignment.GlobalAlphabet._
import modiphy.util.Memo
import scala.collection.LinearSeq


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
    val model = new BasicLikelihoodModel(WAG.pi,WAG.S)

    (0 until 50).foreach{i=>
      val lkl = new SimpleLikelihoodCalc(tree,model,engine=IndexedSeqLikelihoodFactory.apply)
      println(lkl.logLikelihood(aln.columns))
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
  val children:Seq[TreePositionDir]
  def neighbours:Seq[RootedTreePosition] = children
  val ancestralTo:Memo[Node,Boolean] = Memo[Node,Boolean]({n=>
    neighbours.foldLeft(false){(bool,neighbour)=>
      bool || n==neighbour.get || neighbour.ancestralTo(n) //this last one could blow the stack....
    } 
  })
  override def hashCode = children.map{_.hashCode}.foldLeft(get.hashCode + 41){(h1,child)=>
    h1 * 41 + child.hashCode
  }
  override def equals(o:Any)={
    o match {
      case that:TreePositionDir=>false
      case that:RootedTreePosition=>that.get==get && that.children==children
      case _ => false
    }
  }
  val leaves:Seq[Leaf]
}
trait TreePositionDir extends RootedTreePosition{
  val upEdge:Edge
  override def hashCode = super.hashCode + upEdge.hashCode
  override def equals(o:Any)={
    o match {
      case that:TreePositionDir=>that.upEdge==upEdge && that.get==get && that.children==children
      case _ => false
    }
  }
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
    val tree =this
    new TreePositionDir{
      val get = n
      lazy val children = edgeMap(n).filter{e=> ! (e same dir)}.map{e=>traverseDown(e from n right,e)}
      val upEdge = dir
      val leaves = tree.leaves(n,dir)
    }
  }

  def traverseDown(n:Node):RootedTreePosition={
    val tree =this
    new RootedTreePosition{
      val get = n
      lazy val children = edgeMap(n).map{e=>traverseDown(e from n right,e)}
      val leaves = tree.leaves(n)
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
    val oldEdge = edges(i)
    val e = oldEdge
    val newEdge = edges(i).copy(dist=d)
    val myNewEdgeMap = edgeMap.updated(e.left,newEdge::edgeMap(e.left).filter{_!=oldEdge}).updated(e.right,newEdge::edgeMap(e.right).filter{_!=oldEdge})
      copy(newEdges = edges.updated(i,newEdge),newEdgeMap = myNewEdgeMap)
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


class MixtureLikelihoodCalc(priors:Seq[Double],tree:Tree,m:Seq[SingleModel],lkl:Option[Seq[SimpleLikelihoodCalc]]=None){
  import scala.actors.Futures.future
  val lklCalc = lkl.getOrElse{m.map{new SimpleLikelihoodCalc(tree,_)}}
  
  def logLikelihood(patterns:Seq[Pattern])={
    if (Parallel.on){
      patterns.map{pattern=>
          future{
           lklCalc.zip(priors).map{t=> t._1.likelihood(pattern) * t._2}.reduceLeft{_+_}
          }
        }.map{f=>math.log(f())}.foldLeft(0.0D){_+_}
     } else {
      patterns.map{pattern=>
         lklCalc.zip(priors).map{t=> t._1.likelihood(pattern) * t._2}.reduceLeft{_+_}
     }.map{f=>math.log(f)}.foldLeft(0.0D){_+_}
   }
  }
}
class SimpleLikelihoodCalc(tree:Tree,m:SingleModel, var cache:Map[RootedTreePosition,Map[Seq[Letter],PartialLikelihoods]] = Map[RootedTreePosition,Map[Seq[Letter],PartialLikelihoods]](),val engine:LikelihoodEngine=DefaultLikelihoodFactory.apply){
  import SimpleLikelihoodCalc._

import engine.combinePartialLikelihoods
import engine.partialLikelihoodCalc
import engine.finalLikelihood

   def cacheLookup(pos:RootedTreePosition,pattern:Seq[Letter])={
     None
     /*
     cache.get(pos) match {
       case None=>None
       case Some(m2)=>m2 get pattern
     }*/
   }
   def cacheAdd(pos:RootedTreePosition,pattern:Seq[Letter],pl:PartialLikelihoods)={
     cache
     /*
     cache.get(pos) match {
       case None=>cache updated (pos,Map(pattern->pl))
       case Some(m2)=>cache updated (pos,m2 updated (pattern,pl))
     }*/
   }

   
    def partialLikelihoods(treePos:RootedTreePosition,p:Pattern):PartialLikelihoods={
      import scala.actors.Futures.future
      val myPatterns = treePos.leaves.map{leaf=>p(leaf.id.get)}
      cacheLookup(treePos,myPatterns).getOrElse{
        val ans = treePos.get match {
          case n:INode=>
            combinePartialLikelihoods(
              treePos.children.toList.map{tp=>
              val plStart = partialLikelihoods(tp,p)
              partialLikelihoodCalc(plStart,m(tp.upEdge)) 
            }
          )
          case l:Leaf=>leafPartialLikelihoods(p(l.id.get))
        }

        cache = cacheAdd(treePos,myPatterns,ans)
        ans
      }
    }

  def likelihood(p:Pattern,root:RootedTreePosition=tree.traverseDown(tree.defaultRoot)):Likelihood={
    finalLikelihood(partialLikelihoods(root,p),m.pi(root.get))
  }
  def likelihoods(p:Seq[Pattern],root:RootedTreePosition=tree.traverseDown(tree.defaultRoot)):Seq[Double]={
    import scala.actors.Futures._
    if (Parallel.on){ 
      p.grouped(300).toList.map{subList=>future{subList.map{pat=>likelihood(pat,root)}}}.map{_()}.flatten
    }
    else {p.map{pat=>likelihood(pat,root)}}
  }

  def logLikelihood(p:Seq[Pattern],root:RootedTreePosition=tree.traverseDown(tree.defaultRoot)):Double={
   likelihoods(p,root).foldLeft(0.0D){_+math.log(_)}
  }
  val leafPartialLikelihoods:Memo[Letter,PartialLikelihoods]=Memo[Letter,PartialLikelihoods](l=>l match {
    case a if (a.isReal) => List.fill(l.alphabet.length)(0.0).updated(l.id,1.0)
    case a => List.fill(l.alphabet.length)(1.0)
  })


  def factory(t:Tree,m:SingleModel,cache:Cache) = new SimpleLikelihoodCalc(t,m,cache)

  def update(t:Tree)={
    val updatedCache = tree.differences(t).foldLeft(cache){(c2,e)=>
        c2.keys.foldLeft(c2){(map,tp)=> 
         if (tp ancestralTo e.left){ // if I am ancestral to edge e then both nodes are my descendents
           map - tp
         }else {
           map
         }
      }
    }
    factory(t,m,updatedCache)
  }

  def update(m:SingleModel)={
    factory(tree,m,Map[RootedTreePosition,Map[Seq[Letter],PartialLikelihoods]]())
  }
}

trait LikelihoodEngine{
  def combinePartialLikelihoods(intermediates:List[PartialLikelihoods]):PartialLikelihoods
  def partialLikelihoodCalc(end:PartialLikelihoods,matrix:Matrix):PartialLikelihoods    
  def finalLikelihood(partial:PartialLikelihoods,pi:IndexedSeq[Double]):Likelihood
}
class ColtLikelihoodCalc extends LikelihoodEngine{
  import SimpleLikelihoodCalc.Cache
  import cern.colt.matrix._
  val func = new cern.colt.function.DoubleDoubleFunction{def apply(x:Double,y:Double)=x*y}

  import modiphy.math.EnhancedMatrix._

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

class IndexedSeqLikelihoodCalc extends LikelihoodEngine{
  import SimpleLikelihoodCalc.Cache
  import cern.colt.matrix._

  import modiphy.math.EnhancedMatrix._

  def combinePartialLikelihoods(intermediates:List[PartialLikelihoods]):PartialLikelihoods = {
    combinePartialLikelihoods(intermediates,List[Double]())
  }
  def combinePartialLikelihoods(intermediates:List[PartialLikelihoods],ans:List[Double]):PartialLikelihoods={
    if (intermediates.head.isEmpty){
      ans.reverse.toList
    }else {
      combinePartialLikelihoods(intermediates.map{_.tail},intermediates.map{_.head}.product::ans)
    }
  }
  def partialLikelihoodCalc(end:PartialLikelihoods,matrix:Matrix):PartialLikelihoods={
    matrix.map{end.dotProduct}.toList
  }
  def finalLikelihood(partial:PartialLikelihoods,pi:IndexedSeq[Double]):Likelihood={
    coltLikelihoods(List(Seq2Vec(partial)),pi).head
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
