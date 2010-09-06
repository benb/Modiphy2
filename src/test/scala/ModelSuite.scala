import org.scalatest.FunSuite
import org.scalatest.matchers.ShouldMatchers._
import modiphy.alignment._
import modiphy.model._
import modiphy.tree._
import modiphy.math.constants._
import modiphy.opt._
import ModelData._

class ModelSuite extends FunSuite {


   
  test("WAG"){
    WAG.pi.length should be (20)
    WAG.S.length should be (20)
    WAG.S.head.length should be (20)
  }

  test("log likelihood of basic model should match PAML") {
    val model = BasicLikelihoodModel(WAG.pi,WAG.S)
    val tree = Tree(treeStr)
    val aln = Fasta(alnStr).parseWith(AminoAcid)
    /*
    (0 to 50).foreach{i=>
      val tree = Tree(treeStr)
      val lkl = new SimpleLikelihoodCalc(tree,model) 
    
      lkl.logLikelihood(aln.columns) should be (-6057.892394 plusOrMinus 0.001)//from PAML
    }*/
    for (factory <- ColtLikelihoodFactory::IndexedSeqLikelihoodFactory::Nil){
      val lkl = new SimpleLikelihoodCalc(tree,model,aln,engine=factory.apply) 
      lkl.logLikelihood should be (-6057.892394 plusOrMinus 0.001)//from PAML
      val tree2 = tree setBranchLength (2,2.0)
      lkl.updated(tree2).logLikelihood should be < (-6057.892394)
      val tree3 = tree setBranchLength (2, tree.getBranchLength(2))
      lkl.updated(tree2).updated(tree3).logLikelihood should be (-6057.892394 plusOrMinus 0.001) 
    }
    
   
      }
  test("Gamma mixture model should match PAML"){
    val model = GammaModel(WAG.pi,WAG.S,0.5,4)
    val tree = Tree(treeStr)
    val aln = Fasta(alnStr).parseWith(AminoAcid)


    val plusF=aln.frequencies // Vector(0.038195,0.070238,0.054858,0.072802,0.037939,0.046398,0.080749,0.048962,0.017175,0.043066,0.085106,0.069726,0.015124,0.046142,0.028198,0.073571,0.044604,0.024096,0.049474,0.053576)
    val modelF = model updatedVec (Pi,plusF,None)//GammaModel(plusF,WAG.S,0.5,4)

    val lkl = new MixtureLikelihoodCalc(Vector.fill(4)(0.25),tree,aln,model)
    lkl.logLikelihood should be (-5808.929978 plusOrMinus 0.001)
    val lkl2 = new MixtureLikelihoodCalc(Vector.fill(4)(0.25),tree,aln,modelF)
    lkl2.logLikelihood should be (-5810.399586 plusOrMinus 0.001)
  }

  test("Gamma mixture model should match PAML (Controlled by OptModel)"){
    val model = GammaModel(WAG.pi,WAG.S,0.5,4)
    val tree = Tree(treeStr)
    val aln = Fasta(alnStr).parseWith(AminoAcid)
    val lkl = new MixtureLikelihoodCalc(Vector.fill(4)(0.25),tree,aln,model)
    val optModel = new OptModel(model,lkl,tree,aln)
    optModel.logLikelihood should be (-5808.929978 plusOrMinus 0.001)
    optModel(Pi)=aln.frequencies
    optModel.logLikelihood should be (-5810.399586 plusOrMinus 0.001)
  }
}


