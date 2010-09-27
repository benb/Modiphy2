package modiphy.tree
import scala.collection.IndexedSeq
sealed trait Node{
  def allMySubTrees:Set[NonRoot]
  def numNodes:Int
  def numLeaves:Int
}
sealed trait NonRoot extends Node{
  def reRoot(parent:Root):Option[Root]
  def allRootedTrees(parent:Root):Set[Root]
  def toRoot(n:NonRoot):Option[Root]
  def branchLengthString(bl:NonRoot=>Double):String
  def drop(sub:NonRoot):NonRoot
}
class Split(val left:NonRoot,val right:NonRoot) extends Ordered[Split]{
  override def equals(other:Any)= other match {
    case Split(`left`,`right`)=>true
    case Split(`right`,`left`)=>true
    case _ => false
  }
  override def toString = "Split( " + left + " , " + right + ")"
  def compare(that:Split)={
    val a = if (left.toString < right.toString){left.toString}else{right.toString}
    a.compare{
      if (that.left.toString < that.right.toString){that.left.toString}else{that.right.toString}
    }
  }
  override def hashCode = (left.hashCode + right.hashCode) * 43
  /*
   Is the other split compatible with this, i.e. does it represent the same 
   split on a tree with either additional or dropped nodes on one side
  */
  def compatible(other:Split)={
    other.left==left || other.right==left || other.left==right || other.right==right
  }
}
object Split{
  def apply(left:NonRoot,right:NonRoot)=new Split(left,right)
  def unapply(s:Split) = Some(s.left,s.right)
}

object SubTree{
  def apply(c1:NonRoot,c2:NonRoot):SubTree=SubTree(Set(c1,c2))
  def apply(c1:NonRoot,c2:NonRoot,c3:NonRoot):SubTree=SubTree(Set(c1,c2,c3))
}
case class SubTree(children:Set[NonRoot]) extends NonRoot{
  def subTrees=children.filter{_ match {case t:SubTree=>true; case _=>false}}.map{_.asInstanceOf[SubTree]}
  def reRoot(parent:Root) = parent.childReRoot(this)
    override def toString="(" + children.mkString(",") + ")"
  def allRootedTrees(oldRoot:Root):Set[Root] = {
    val meRoot = oldRoot.childReRoot(this).get
    val local = children.map{_.reRoot(meRoot)}.filterNot{_==None}.map{_.get}
    children.foldLeft(local){(m,c)=> m ++ c.allRootedTrees(meRoot)}
  }
  //FIXME
  def toRoot(c3:NonRoot)=Some(Root(children + c3))


  def allMySubTrees = children.map{_.allMySubTrees}.foldLeft(children){(m,a) => m ++ a}// + this
  def branchLengthString(bl:NonRoot=>Double)={
    "(" + children.map{_.branchLengthString(bl)}.mkString(",") + "):" + bl(this).toString
  }
  def drop(sub:NonRoot)={
    val newChildren = children.filterNot{_==sub}.map{_.drop(sub)}.map{node => 
     node match {
       case l:Leaf => l
       case s:SubTree => 
         if (s.children.size==1){
           s.children.head
         }else {
           node
         }
       }
     }
     SubTree(newChildren)
   }
  def numNodes = 1 + children.foldLeft(0){_+_.numNodes}
  def numLeaves = children.foldLeft(0){_+_.numLeaves}

}
object Root{
  def apply(c1:NonRoot,c2:NonRoot,c3:NonRoot):Root=Root(Set(c1,c2,c3))
}
case class Root(children:Set[NonRoot]) extends Node{
  def childReRoot(n:NonRoot):Option[Root] = n match {
    case t:SubTree if (children(t))=> t.toRoot(SubTree(children-t))
    case _ => None
  }
  def allRootedChildren:Set[Root] = children.map{childReRoot}.filterNot{_==None}.map{_.get}
  lazy val allRootedTrees:Set[Root] = {
    val start:Set[Set[Root]] = children.map{_.allRootedTrees(this)}
    val ans = start.foldLeft(allRootedChildren + this){_++_}
    assert(ans.size > 1)
    ans
  }

  def allMySubTrees:Set[NonRoot] = children.flatMap{_.allMySubTrees} ++ children
  /*
   Seq of all possible SubTrees from all possible rootings
  */
  lazy val allSubTrees:IndexedSeq[NonRoot]={
    val ans = allRootedTrees.map{_.allMySubTrees}.foldLeft(Set[NonRoot]()){_++_}.toIndexedSeq
    assert(ans.length > 1)
    ans
  }
  lazy val allSplits:IndexedSeq[Split]=
   allRootedTrees.toList.flatMap{root=> root.children.map{x=>
     if (root.children.size >2){
      Split(x,SubTree(root.children-x))
    }else {
      Split(x,(root.children-x).head)
    }
    }
   }.toIndexedSeq.sorted.distinct
  lazy val split:Map[NonRoot,Split]={
    allSplits.flatMap{s => (s.left,s)::(s.right,s)::Nil}.toMap
  }
  def reRoot(i:Int):Root=reRoot(allSplits(i).left)
  def reRoot(anySub:NonRoot):Root={ allRootedTrees.find{root=> root.children.contains(anySub)}.get }
  override def toString="(" + children.mkString(",") + ");"
  def branchLengthString(blSeq:IndexedSeq[Double]):String={
    branchLengthString(allSplits.map{_.left}.zip(blSeq).toMap ++ allSplits.map{_.right}.zip(blSeq))
  }
  def branchLengthString(bl:Map[NonRoot,Double]):String={
    "(" + children.map{_.branchLengthString(bl)}.mkString(",") + ");"
  }
  def drop(sub:NonRoot)={
    val newChildren = children.filterNot{_==sub}.map{_.drop(sub)}.map{node => 
     node match {
       case l:Leaf => l
       case s:SubTree => 
         if (s.children.size==1){
           s.children.head
         }else {
           s
         }
       }
     }
     Root(newChildren)
  }
  def numNodes = 1 + children.foldLeft(0){_+_.numNodes}
  def numLeaves = children.foldLeft(0){_+_.numLeaves}
}
case class Leaf(name:String) extends NonRoot{
  override def toString = name
  def reRoot(n:Root)=None
  def toRoot(n:NonRoot)=None
  def allRootedTrees(parent:Root)=Set()
  def allMySubTrees = Set()
  def branchLengthString(bl:NonRoot=>Double)={
    name + ":" +  bl(this).toString
  }
  def drop(sub:NonRoot)=this
  def numNodes=1
  def numLeaves=1
}

object TreeTest{
  implicit def toLeaf(s:String)=Leaf(s)
  implicit def toSubTree(t:(NonRoot,NonRoot))=SubTree(t._1,t._1)
  def main(args:Array[String]){
    val t1 = Root(SubTree(SubTree("onea","oneb"),"two"),SubTree(SubTree("threea","threeb"),"four"),SubTree("five","six"))
    println(t1)
    println(t1.allRootedChildren)
    println(t1.allSubTrees)
    println(t1.allSplits.map{t=> t.left.toString + "/" + t.right.toString})
    println(t1.allRootedTrees.mkString("\n"))
    println(t1.branchLengthString((0.0 to 2.0 by 0.1).toIndexedSeq))
    println(t1.allRootedChildren.toList(2).branchLengthString((0.0 to 2.0 by 0.1).toIndexedSeq))
  }
}

case class Tree(root:Root,bl:IndexedSeq[Double]){
  def getBranchLengths = bl
  def setBranchLengths(bl2:IndexedSeq[Double])=copy(root,bl2)
  def setBranchLength(i:Int,bl2:Double)=copy(root,bl.updated(i,bl2))
  override def toString = root.branchLengthString(bl)
  lazy val treeLength = bl.reduceLeft(_+_)
  def drop(sub:NonRoot):Tree={
    val newRoot = root.drop(sub)
    val newBl = root.allSplits.map{s1=>newRoot.allSplits.find(s2=> s2 compatible s1)}.zip(bl).foldLeft(Map[Split,Double]()){(m,t)=> val (newSplit,bl)=t
      if (newSplit.isDefined){
        m.updated(newSplit.get,m.getOrElse(newSplit.get,0.0)+bl) 
      }else {
        m
      }
    }.toIndexedSeq
    val newBl2 = newBl.sortBy{_._1}.map{_._2}
    Tree(newRoot,newBl2)
  }
  def drop(s:String):Tree={
   // val sub = new TreeParser{def parseAll=parse(node,s)}.parseAll.get._1
  //  drop(sub)
    drop(Leaf(s))
  }
  lazy val branchLength:Map[NonRoot,Double]={
    root.allMySubTrees.map{root.split}.toList.sorted.zip(bl).flatMap{t=>List((t._1.left,t._2),(t._1.right,t._2))}.toMap
  }
  def branchLength(i:Int):Double=bl(i)
  def numNodes = root.numNodes
  def numLeaves = root.numLeaves
  def makeTrifurcating={
    assert(root.children.size==2)
    val left:SubTree = root.children.find{_.isInstanceOf[SubTree]}.get.asInstanceOf[SubTree]
    val right:NonRoot = (root.children-left).head
    val newRootBL = branchLength(left) + branchLength(right)
    val newRoot = Root(left.children + right)
    val blMap = root.allSplits.zip(bl).toMap 
    val newBLMap = (blMap.updated(root.split(right),blMap(root.split(right))+blMap(root.split(left))) - root.split(left)).toIndexedSeq

    println(blMap)
    println(newBLMap)

    Tree(newRoot, newBLMap.sortBy(_._1).map{_._2})
  }
}
object Tree{
  def apply(newick:String):Tree = {
    new TreeParser{def parseAll=parse(tree,newick)}.parseAll.get
  }
}
