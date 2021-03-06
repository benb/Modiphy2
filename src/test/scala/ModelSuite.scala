import org.scalatest.FunSuite
import org.scalatest.matchers.ShouldMatchers._
import modiphy.alignment._
import modiphy.model._
import modiphy.tree._
import modiphy.math.constants._
import modiphy.opt._
import modiphy.calc._
import modiphy.math.EnhancedMatrix._
import ModelData._
import modiphy.model.Types._

class ModelSuite extends FunSuite {


   
  test("WAG"){
    WAG.pi.length should be (20)
    WAG.S.length should be (20)
    WAG.S.head.length should be (20)
  }

  test("appliesTo"){
    val model = new BasicLikelihoodModel(List(Pi << WAG.pi,S << WAG.S,Rate << 1.0),2)
    model.appliesToMe(Pi,MatchP(2)) should be (true)
    model.appliesToMe(Pi,MatchP(1)) should be (false)
    model.appliesToMe(Pi,MatchSet(0,1)) should be (false)
    model.appliesToMe(Pi,MatchSet(2,3)) should be (true)
    model.appliesToMe(S,MatchAll) should be (true)
    model.appliesToMe(Gamma,MatchAll) should be (false)
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
      val tree3 = tree setBranchLength (2, tree.getBranchLengths(2))
      lkl.updated(tree2).updated(tree3).logLikelihood should be (-6057.892394 plusOrMinus 0.001) 
    }
  }
  test("Gamma mixture model should match PAML"){
    val model = GammaModel(WAG.pi,WAG.S,0.5,4)
    val tree = Tree(treeStr)
    val aln = Fasta(alnStr).parseWith(AminoAcid)


    val plusF=aln.frequencies // Vector(0.038195,0.070238,0.054858,0.072802,0.037939,0.046398,0.080749,0.048962,0.017175,0.043066,0.085106,0.069726,0.015124,0.046142,0.028198,0.073571,0.044604,0.024096,0.049474,0.053576)
    val modelF = model updatedVec (Pi,plusF,MatchAll)//GammaModel(plusF,WAG.S,0.5,4)

    val lkl = new MixtureLikelihoodCalc(tree,aln,model)
    lkl.logLikelihood should be (-5808.929978 plusOrMinus 0.001)
    val cmp = lkl.componentLikelihoods.asInstanceOf[List[List[Double]]]// should work
    cmp.length should be (4)
    cmp.reduceLeft{(a,b)=>
      a.length should be (b.length)
      a
    }

    val flt = lkl.flatComponentLikelihoods
    flt.length should be (4)
    flt.reduceLeft{(a,b)=>
      a.length should be (b.length)
      a
    }
    println(lkl.flatComponentLikelihoods.map{_.length})
    val lkl2 = new MixtureLikelihoodCalc(tree,aln,modelF)
    lkl2.logLikelihood should be (-5810.399586 plusOrMinus 0.001)
  }

  test("Gamma mixture model should match PAML (Controlled by OptModel)"){
    val model = GammaModel(WAG.pi,WAG.S,0.5,4)
    val tree = Tree(treeStr)
    val aln = Fasta(alnStr).parseWith(AminoAcid)
    val lkl = new MixtureLikelihoodCalc(tree,aln,model)
    val optModel = new OptModel(lkl,tree,aln)
    optModel.logLikelihood should be (-5808.929978 plusOrMinus 0.001)
    optModel(Pi)=aln.frequencies
    optModel.logLikelihood should be (-5810.399586 plusOrMinus 0.001)
  }

  test("Opt Gamma Model should match PAML"){
    val tree = Tree(treeStr)
    val aln = Fasta(alnStr).parseWith(AminoAcid)
    val model = GammaModel(aln.frequencies,WAG.S,0.5,4)
    val lkl = new MixtureLikelihoodCalc(tree,aln,model)
    val optModel = new OptModel(lkl,tree,aln)
    val t1 = System.currentTimeMillis
    optModel optimiseAll Gamma
    val t2 = System.currentTimeMillis
    println(t2-t1)
    optModel.logLikelihood should be > ( -5809.180030)
    println(optModel.logLikelihood)
    println(optModel(Gamma))
    optModel(Gamma).get.head should be (0.57932 plusOrMinus 0.01)
   // optModel optimiseAll BranchLengths
  }
  test{"WTF"}{
    val tree = Tree(treeStr)
    val aln = Fasta(alnStr).parseWith(AminoAcid)
    val zeroModel = BasicLikelihoodModel(WAG.pi,WAG.S,3.0,1)
   // val modelG = GammaModel(WAG.pi,WAG.S,0.5,4) add (zeroModel,0.99)
    val modelG = zeroModel.updatedRate(3.0) add (zeroModel,0.5,1) updatedRate 3.0
    println("Fail " + modelG.rate)
    println("Fail " + zeroModel.rate)
  
    val lkl = new MixtureLikelihoodCalc(tree,aln,modelG)
    val lkl2 = new SimpleLikelihoodCalc(tree,zeroModel,aln)
    val lkl3 = new SimpleLikelihoodCalc(tree,modelG,aln)
    lkl2.logLikelihood should be (lkl.logLikelihood plusOrMinus 0.0001)
    lkl3.logLikelihood should be (lkl.logLikelihood plusOrMinus 0.001)
//    lkl2.logLikelihood should be (-5864.879865 plusOrMinus 0.01) // phyml -d aa -o n -i png1-aln.phy -u png1.tre -a 0.5 -m WAG -v 0.2
//  lkl.logLikelihood should be (-5864.879865 plusOrMinus 0.01) 
    
  }
  
  test{"Gamma+I"}{
    val tree = Tree(treeStr)
    val aln = Fasta(alnStr).parseWith(AminoAcid)
    val zeroModel = BasicLikelihoodModel.zeroRate(WAG.pi,1)
    val modelG = GammaModel(WAG.pi,WAG.S,0.5,4) add (zeroModel,0.4)
    println("Fail " + modelG.rate)
    println("Fail " + zeroModel.rate)
  
    val lkl = new MixtureLikelihoodCalc(tree,aln,modelG)
    val opt1 = new OptModel(lkl,tree,aln)
    val lkl2 = new SimpleLikelihoodCalc(tree,modelG,aln)
    val opt2 = new OptModel(lkl2,tree,aln)
    lkl2.logLikelihood should be (lkl.logLikelihood plusOrMinus 0.0001)
    lkl2.logLikelihood should be (-5864.879865 plusOrMinus 0.01)
    opt2(MixturePrior)=Vector(0.8,0.2)
    opt1(MixturePrior)=Vector(0.8,0.2)
    opt2.logLikelihood should be (-5824.746968 plusOrMinus 0.01)
    opt1.logLikelihood should be (-5824.746968 plusOrMinus 0.01)

//    lkl2.logLikelihood should be (-5864.879865 plusOrMinus 0.01) // phyml -d aa -o n -i png1-aln.phy -u png1.tre -a 0.5 -m WAG -v 0.2
//  lkl.logLikelihood should be (-5864.879865 plusOrMinus 0.01) 
    
  }
  test("Site class model with 1 class"){
    val tree = Tree(treeStr)
    val aln = Fasta(alnStr).parseWith(AminoAcid)
    val model = BasicLikelihoodModel(WAG.pi,WAG.S)
    val thmm = StdSiteClassModel(List(model,model))
    val lkl = new SimpleLikelihoodCalc(tree,thmm,aln)
    val optModel = new OptModel(lkl,tree,aln)
    optModel.logLikelihood should be (-6057.892394 plusOrMinus 0.001)
  }


  test("Opt Gamma Model should match PAML (thmm)"){
    val tree = Tree(treeStr)
    val aln = Fasta(alnStr).parseWith(AminoAcid)
    val model = GammaModel(aln.frequencies,WAG.S,0.5,4)
    val lkl = new SimpleLikelihoodCalc(tree,model,aln)
    val optModel = new OptModel(lkl,tree,aln)
    optModel.logLikelihood should be (-5810.399586 plusOrMinus 0.001)
  }
  test("THMM"){
    val tree = Tree(pfTree)
    val aln = Fasta(pfAln).parseWith(AminoAcid)
    val pi = aln.frequencies//Vector(0.024191,0.002492,0.002932,0.002492,0.001906,0.002492,0.006304,0.023018,0.002346,0.026683,0.034307,0.008943,0.007037,0.014808,0.005278,0.018326,0.013928,0.007477,0.007917,0.020379).normalise
    println("PI " + pi)
    val zeroModel = BasicLikelihoodModel.zeroRate(pi,1)
    val gammaModel = GammaModel(pi,WAG.S,3.270690,4) add (zeroModel,0.066963)
    val thmm = new ThmmSiteClassModel((Sigma << 2.415327)::Nil, 1,gammaModel,AminoAcid.matLength )
    val lkl = new SimpleLikelihoodCalc(tree,thmm,aln)
    val optModel = new OptModel(lkl,tree,aln)
    println(thmm.params)
    optModel.logLikelihood should be (-2973.2188766283607 plusOrMinus 1e-2) // lnL from modiphy1
    //optModel.optimise((MixturePrior,Some(1)))
  }
}


