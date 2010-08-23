package modiphy.tree

import scala.util.parsing.combinator._
/**
  Used internally for parsing trees
*/

class TreeParser() extends JavaTokenParsers{
  type BranchLength=(Double,Option[String])

  def tree: Parser[Tree] = "("~node~rep(","~>node)~");" ^^ {
    case "("~node1~nodeList~");" => Tree({val myNode = new INode; (node1::nodeList).map{t=> Edge(myNode,t._1,t._2._1)} ++ (node1::nodeList).map{_._3}.flatten})
  }

  def node: Parser[(Node,BranchLength,List[Edge])] = leaf | "("~node~rep(","~>node)~"):"~branchLength ^^ {
    case "("~node1~nodeList~"):"~length => {val myNode = new INode; (myNode,length,(node1::nodeList).map{t=> Edge(myNode,t._1,t._2._1)} ++ (node1::nodeList).map{_._3}.flatten)}
  }
  
  def seqName: Parser[String] = regex(new scala.util.matching.Regex("[a-zA-Z0-9_.+-]+"))

  def leaf: Parser[(Leaf,BranchLength,List[Edge])] = seqName~":"~branchLength ^^ {case name~":"~branchLength => (Leaf(name),branchLength,Nil)}

  def branchLength:Parser[BranchLength] = {floatingPointNumber~"#"~seqName ^^{
    case length~"#"~id => (length.toDouble match {case c if c<0.0 => 0.0; case c=>c},Some(id))
  } | floatingPointNumber  ^^ {
    case length => (length.toDouble match {case c if c<0.0 => 0.0; case c => c},None)
  }
  }
}
