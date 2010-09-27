package modiphy.tree

import scala.util.parsing.combinator._
import scala.collection.immutable._
/**
  Used internally for parsing trees
*/

class TreeParser() extends JavaTokenParsers{
  type BranchLength=(Double,Option[String])
  type BranchLengthMap = Map[NonRoot,BranchLength]
  type BranchLengthFinalMap = Map[Split,BranchLength]

  def tree: Parser[Tree] = "("~node~rep(","~>node)~");" ^^ {
    case "("~node1~nodeList~");" => 
      val fullList:List[(NonRoot,BranchLengthMap)] = node1::nodeList
      val root = Root(fullList.map{_._1}.toSet)
      val finBL:IndexedSeq[(Split,BranchLength)] = fullList.map{_._2}.reduceLeft(_++_).map{t => (root.split(t._1),t._2)}.toIndexedSeq
      val finBL2 = finBL.sortBy(_._1).map{_._2._1}
      assert(finBL2.length == root.allMySubTrees.size)
      Tree(root,finBL2)
  }

  def node: Parser[(NonRoot,BranchLengthMap)] = leaf | "("~node~rep(","~>node)~"):"~branchLength ^^ {
    case "("~node1~nodeList~"):"~length => {
      val fullList = node1::nodeList
      val nonRoot = SubTree(fullList.map{_._1}.toSet)
      val blMap = fullList.map{_._2}.reduceLeft(_++_) + ((nonRoot,length))
      (nonRoot,blMap)
    }
  }
  
  def seqName: Parser[String] = regex(new scala.util.matching.Regex("[a-zA-Z0-9_.+-]+"))

  def leaf: Parser[(NonRoot,BranchLengthMap)] = seqName~":"~branchLength ^^ {case name~":"~branchLength => (Leaf(name),Map(Leaf(name)->branchLength))}

  def branchLength:Parser[BranchLength] = {floatingPointNumber~"#"~seqName ^^{
    case length~"#"~id => (length.toDouble match {case c if c<0.0 => 0.0; case c=>c},Some(id))
  } | floatingPointNumber  ^^ {
    case length => (length.toDouble match {case c if c<0.0 => 0.0; case c => c},None)
  }
  }
}
