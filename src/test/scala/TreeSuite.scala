import org.scalatest.FunSuite
import org.scalatest.matchers.ShouldMatchers._
import modiphy.tree._

class TreeSuite extends FunSuite {
  val newick = "(Chlamydomona:0.082441#0,Chlorella:0.050529#1,((Euglena:0.099344#5,Nephroselmis:0.088455#10):0.010341#27,(((Mycobacteriu:0.527327#8,(Mycoplasma:0.338658#9,(Ecoli:0.214706#4,Rickettsia:0.253023#13):0.050486#16):0.061350#17):0.143567#18,(Synechocysti:0.079101#14,(Thermosynech:0.057044#15,(Gloeobacter:0.179811#6,Prochlorococ:0.109110#12):0.038575#19):0.009400#20):0.043594#21):0.033610#22,(Guillardia:0.067446#7,(Porphyra:0.064926#11,(Cyanidioschy:0.076726#2,Cyanidium:0.107181#3):0.028576#23):0.013723#24):0.010431#25):0.060348#28):0.000000#26);"
  val tree = Tree(newick)

  test("parse"){
    println(tree)
    val otherNode = tree.traverseFrom("Euglena").get.neighbours.map{_.get}.collect{case i:INode=>i}.head
    tree.getBranchLengths.length should equal (29)

    println(tree reRoot otherNode)
    val bl1 = tree branchLength 3
    val t2 = tree setBranchLength(3,bl1+1.3)
    t2.treeLength should equal (tree.treeLength+1.3)
    println(t2)
  }
  test("Edit"){
    val tree2 = tree.drop("Euglena")
    tree2 should not be(tree)
    tree2.numNodes should be (tree.numNodes - 2)
    tree2.numLeaves should be (tree.numLeaves - 1)
    println(tree2)
  }
}
