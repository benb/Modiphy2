import org.scalatest.FunSuite
import org.scalatest.matchers.ShouldMatchers._
import modiphy.tree._

class TreeSuite extends FunSuite {
//  val newick = "(((One:1,Two:2):3,Nine:9):12,(Four:4,Five:5):9,Six:6);"
  val newick = "(Chlamydomona:0.082441#0,Chlorella:0.050529#1,((Euglena:0.099344#5,Nephroselmis:0.088455#10):0.010341#27,(((Mycobacteriu:0.527327#8,(Mycoplasma:0.338658#9,(Ecoli:0.214706#4,Rickettsia:0.253023#13):0.050486#16):0.061350#17):0.143567#18,(Synechocysti:0.079101#14,(Thermosynech:0.057044#15,(Gloeobacter:0.179811#6,Prochlorococ:0.109110#12):0.038575#19):0.009400#20):0.043594#21):0.033610#22,(Guillardia:0.067446#7,(Porphyra:0.064926#11,(Cyanidioschy:0.076726#2,Cyanidium:0.107181#3):0.028576#23):0.013723#24):0.010431#25):0.060348#28):0.000000#26);"
  val tree = Tree(newick)

  test("Split"){
    val split = tree.root.allSplits.head
    split should equal (Split(split.right,split.left))
  }
  test("parse1"){
    println(tree)
    tree.getBranchLengths.length should equal (29)
    val bl1 = tree branchLength 3
    val t2 = tree setBranchLength(3,bl1+1.3)
    t2.treeLength should equal (tree.treeLength+1.3)
    val t3 = t2.setBranchLengths(Vector.fill(29){0.1})
    t3.getBranchLengths.length should equal (tree.getBranchLengths.length)
    println(t2)
    println(t3)
  }
  test("Nodes"){
    val tree = Tree("(((One:1,Two:2):3,Nine:9):12,(Four:4,Five:5):9,Six:6);")
    tree.numLeaves should be (6)
    tree.numNodes should be (10)
    tree.branchLength(Leaf("Six")) should be (6.0)
  }
  test("Nodes2"){
    val tree = Tree("(((One:1,Two:2):3,Nine:9):10,((Four:4,Five:5):9,Six:6):2);")
    println(tree)
    tree.numLeaves should be (6)
    tree.numNodes should be (10)
    tree.branchLength(Leaf("Six")) should be (6.0)
  }

  test("Edit1"){
    val tree2 = tree.drop("Euglena")
    tree2 should not be(tree)
    tree2.numNodes should be (tree.numNodes - 2)
    tree2.numLeaves should be (tree.numLeaves - 1)
    println(tree2)
  }
  test("Edit2"){
    val tree2 = Tree("(((((((3:0.00553833,6:0.00585932):0.00091763,7:0.00638805):0.00073745,(8:0.00654679,9:0.00672181):0.00103863):0.01523047,(12:0.00309931,(11:0.00385310,(1:0.00348553,10:0.00406893):0.00006191):0.00011043):0.01185213):0.00594087,23:0.00799371):0.0141801,(5:0.00641784,(19:0.00701802,((21:0.00000801,24:0.00008238):0.00239273,(15:0.00209236,(2:0.00085555,17:0.00054690):0.00116256):0.00049543):0.00097509):0.00197679):0.00320602):0.01813072,(25:0.00492161,(18:-0.00000185,20:0.00004704):0.00919822):0.00474103,(14:0.01846822,(28:0.00634597,(27:0.00728544,((4:0.00323671,26:0.00333289):0.00497964,(13:0.00542113,(30:0.00519839,(22:0.00624519,29:0.00556161):0.00061215):0.00182289):0.00103519):0.00071540):0.00105383):0.00613338):0.07993945);")
    val tree3 = tree2.restrictTo(List(12, 11, 10, 1, 9, 8, 7, 6, 3, 25, 24, 23, 21, 20, 19, 18, 17, 15, 5, 2).map{_.toString})
    val tree3a = Tree("(((((((3:0.00553833,6:0.00585932):0.00091763,7:0.00638805):0.00073745,(8:0.00654679,9:0.00672181):0.00103863):0.01523047,(12:0.00309931,(11:0.00385310,(1:0.00348553,10:0.00406893):0.00006191):0.00011043):0.01185213):0.00594087,23:0.00799371):0.0141801,(5:0.00641784,(19:0.00701802,((21:0.00000801,24:0.00008238):0.00239273,(15:0.00209236,(2:0.00085555,17:0.00054690):0.00116256):0.00049543):0.00097509):0.00197679):0.00320602):0.01813072,(25:0.00492161,(18:-0.00000185,20:0.00004704):0.00919822):0.00474103);")
    tree3a.numLeaves should be (tree3.numLeaves)
    tree3a.numNodes should be (tree3.numNodes)
    tree3a.root.allRootedTrees should be (tree3.root.allRootedTrees)
    tree3a.root.allSplits should be (tree3.root.allSplits)
    println(tree3a.bl)
    println(tree3.bl)
    tree3a.bl.zip(tree3.bl).foreach{ t=>
      t._1 should be (t._2 plusOrMinus 0.0001)
    }

    tree2.restrictTo(List(9, 8, 7, 6, 3).map{_.toString})

  }
}
