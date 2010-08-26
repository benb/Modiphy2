import org.scalatest.FunSuite
import org.scalatest.matchers.ShouldMatchers._
import modiphy.alignment._
import modiphy.model._
import modiphy.tree._
import modiphy.math.constants._
import ModelData._

class ModelSuite extends FunSuite {


   
  test("WAG"){
    WAG.pi.length should be (20)
    WAG.S.length should be (20)
    WAG.S.head.length should be (20)
  }

  test("log likelihood of basic model should match PAML") {
    val model = new BasicLikelihoodModel(WAG.pi,WAG.S)
    val tree = Tree(treeStr)
    val aln = Fasta(alnStr).parseWith(AminoAcid)
    val lkl = new SimpleLikelihoodCalc(tree,model) with ColtLikelihoodCalc
    
    lkl.logLikelihood(aln.columns) should be (-6057.892394 plusOrMinus 0.001)//from PAML
    val tree2 = tree setBranchLength (2,2.0)
    lkl update tree2 logLikelihood aln.columns should be < (-6057.892394)
    val tree3 = tree setBranchLength (2, tree.getBranchLength(2))
    println(lkl update tree2 logLikelihood aln.columns)
    lkl update tree2 update tree3 logLikelihood aln.columns should be (-6057.892394 plusOrMinus 0.001) 
  }


}


