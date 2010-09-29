package modiphy.calc
import modiphy.tree._
import scala.collection.immutable.Map.WithDefault
import modiphy.alignment._
import modiphy.model._
import modiphy.alignment.GlobalAlphabet._
import scala.collection.LinearSeq

import modiphy.opt._

object Parallel{
  var on = true
  import jsr166y._
  var forkJoinPool = new ForkJoinPool
  lazy val threshold = -1
  private[this] var nt=0
  def numThreads_=(nT:Int){
    nt=nT
    forkJoinPool = new ForkJoinPool(nT)
  }
  def numThreads=nt
}

object LikelihoodTypes{
  type Pattern=Int=>Letter
  type PartialLikelihoods = LinearSeq[Double]
  type Likelihood = Double
  type LogLikelihood = Double
}
import LikelihoodTypes._

object SimpleLikelihoodCalc{
  type Cache=Map[Node,Map[Seq[Letter],PartialLikelihoods]]
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

class MixtureLikelihoodCalc(tree:Tree,aln:Alignment,m:Model,lkl:Option[Seq[LikelihoodCalc]]=None) extends LikelihoodCalc{
  def priors = m.priors
  val models = m.models
  val lklCalc = lkl.getOrElse{models.map{_.likelihoodCalc(tree,aln)}}
  
  lazy val subLikelihoods  = {
    val myLikelihoods = if (Parallel.on && aln.patternLength > 1000){
      import jsr166y._
      class Calc(subModels:List[(LikelihoodCalc,Double)]) extends RecursiveTask[Seq[Seq[Double]]]{
        def compute = {
          subModels match {
            case model::Nil => subModels.map{t=> t._1.likelihoods.map{_ * t._2}}
            case model::tail => subModels.map{model => new Calc(List(model))}.map{_.fork}.map{_.join}.map{_.head}
            case Nil => Nil
          }
        }
      }
      val calc = new Calc(lklCalc.zip(priors).toList)
      Parallel.forkJoinPool submit calc
      calc.join
    }else {
      lklCalc.zip(priors).map{t=> 
        t._1.likelihoods.map{_ * t._2}
      }
    }
    myLikelihoods
  }
  /*
    list of likelihoods for each component in order
    each component will supply a type A, where A is either
    Seq[Double] or Seq[A]. If there are many nested submodels then
    the rabbit hole could go quite deep...

    It is intended to use pattern matching to extract the A's
  */
  lazy val componentLikelihoods:List[List[_]] = {
    def scale(prior:Double,seq:List[_]):List[_]={
      seq.map{sub=>sub match {
        case d:Double=>{d*prior}
        case s:List[_] =>scale(prior,s)
      }}
    }
    lklCalc.zip(priors).map{t=> scale(t._2,t._1.componentLikelihoods) }.toList
  }
  lazy val flatComponentLikelihoods:List[List[Double]]={
    def flatten(ans:List[List[Double]],in:List[List[_]]):List[List[Double]]={
      in match {
        case Nil => ans
        case s => s.head.head match {
          case s2:List[_] => flatten(flatten(ans,s.head.asInstanceOf[List[List[_]]]),s.tail)
          case d:Double => flatten(s.head.asInstanceOf[List[Double]]::ans,s.tail)
        }
      }
    }
    flatten(Nil,componentLikelihoods).reverse
  }
  
  lazy val posteriors:Seq[Seq[Double]] = {
    val iters = flatComponentLikelihoods.map{_.iterator}
    var ans = List[List[Double]]()
    while (iters.head.hasNext){
      val p = iters.map{_.next}
      val t = p.reduceLeft{_+_}
      ans = p.map{_/t}::ans
    }
    ans.reverse
  }
  lazy val likelihoods = {
    val myLikelihoods = subLikelihoods.map{_.iterator}
    var ans = List[Double]()
    while (myLikelihoods.head.hasNext) {
      ans = myLikelihoods.map{_.next}.reduceLeft{_+_} ::ans
    }
    ans.reverse
  }
  lazy val logLikelihood={
    val patternCount = aln.countList
    likelihoods.zip(patternCount).map{t=> math.log(t._1)*t._2}.reduceLeft{_+_}
  }

  def updated(t:Tree)=new MixtureLikelihoodCalc(t,aln,m,Some(lklCalc.map{_.updated(t)}))
  def updated(m:Model)={
    val newModelList = m.models
    val newLkl = lklCalc.zip(newModelList).map{t=> t._1 updated t._2}
    new MixtureLikelihoodCalc(tree,aln,m,Some(newLkl))
  }

  def updatedVec(p:VectorParamName,vec:IndexedSeq[Double],paramIndex:ParamMatcher)={ 
    updated(m.updatedVec(p,vec,paramIndex))
  }

  def updatedSingle(p:SingleParamName,d:Double,paramIndex:ParamMatcher)={ updated(m.updatedSingle(p,d,paramIndex)) }
  def updatedMat(p:MatrixParamName,mat:IndexedSeq[IndexedSeq[Double]],paramIndex:ParamMatcher)={ updated(m.updatedMat(p,mat,paramIndex)) }
  def model =m
  def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:ParamMatcher)={ updated(m.setOptParam(p,vec,paramIndex))}
}

trait LikelihoodCalc{
   def updated(t:Tree):LikelihoodCalc
   def updated(m:Model):LikelihoodCalc
   def logLikelihood:Double
   def updatedVec(p:VectorParamName,vec:IndexedSeq[Double],paramIndex:ParamMatcher):LikelihoodCalc
   def updatedMat(p:MatrixParamName,mat:IndexedSeq[IndexedSeq[Double]],paramIndex:ParamMatcher):LikelihoodCalc
   def updatedSingle(p:SingleParamName,d:Double,paramIndex:ParamMatcher):LikelihoodCalc
   def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:ParamMatcher):LikelihoodCalc
   def model:Model
   def likelihoods:LinearSeq[Double]
   def posteriors:Seq[Seq[Double]]
   def componentLikelihoods:List[List[_]]
   def flatComponentLikelihoods:LinearSeq[LinearSeq[Double]]
}
class SimpleLikelihoodCalc(val tree:Tree,m:Model,val aln:Alignment,val engine:LikelihoodEngine=DefaultLikelihoodFactory.apply) extends LikelihoodCalc{
  def this(tree:Tree,aln:Alignment,m:Model)=this(tree,m,aln)
  import SimpleLikelihoodCalc._
  
  import engine.combinePartialLikelihoods
  import engine.partialLikelihoodCalc
  import engine.finalLikelihood

  val numClasses = m.numClasses

  
  def componentLikelihoods=List(likelihoods.toList)
    def posteriors=List(likelihoods.map{f=>1.0})
  def model=m
  def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:ParamMatcher)={ updated(m.setOptParam(p,vec,paramIndex))}

   def updatedVec(p:VectorParamName,vec:IndexedSeq[Double],paramIndex:ParamMatcher)={ updated(m.updatedVec(p,vec,paramIndex)) }
   def updatedSingle(p:SingleParamName,d:Double,paramIndex:ParamMatcher)={ updated(m.updatedSingle(p,d,paramIndex)) }
   def updatedMat(p:MatrixParamName,mat:IndexedSeq[IndexedSeq[Double]],paramIndex:ParamMatcher)={ updated(m.updatedMat(p,mat,paramIndex)) }
  def realPartialLikelihoodCalc(treePos:Node):LinearSeq[PartialLikelihoods]={
    val p = aln.patternList
       treePos match {
          case iNode:SubTree=>
           partialLikelihoodCalc(
            combinePartialLikelihoods(iNode.children.toList.map{realPartialLikelihoodCalc}
           ),m(tree,iNode))
          case root:Root=>
            combinePartialLikelihoods(root.children.toList.map{realPartialLikelihoodCalc})
          case l:Leaf=>   
          partialLikelihoodCalc(
            p.map{p2=>
              leafPartialLikelihoods(p2(aln.seqId(l.name)))},
              m(tree,l)
            )
        }
    }

   def flatComponentLikelihoods:LinearSeq[LinearSeq[Double]] = List(likelihoods) 


    import scalaz.Scalaz._

    def partialLikelihoods(treePos:Node)={
      Parallel.on match {
        case true => parallelPartialLikelihoods(treePos)
        case false => realPartialLikelihoodCalc(treePos)
      }
    }
    def parallelPartialLikelihoods(treePos:Node)={
      import jsr166y._
      import Parallel._
      val p = aln.patternList
      class Calc(treePos:Node) extends RecursiveTask[LinearSeq[PartialLikelihoods]]{
        def compute():LinearSeq[PartialLikelihoods]={
          treePos match {
            case sub:SubTree=>
                val childCalc = sub.children.toList.map{new Calc(_)}.map{_.fork}
                partialLikelihoodCalc( combinePartialLikelihoods(childCalc.map{_.join}) , m(tree,sub) )
            case root:Root=>
              val childCalc = root.children.toList.map{new Calc(_)}.map{_.fork}
              combinePartialLikelihoods( childCalc.map{_.join})
            case l:Leaf=> partialLikelihoodCalc(allLeafPartialLikelihoods(p,l),m(tree,l))
          }
        }
      }
      val calc = new Calc(treePos)
      forkJoinPool submit calc
      calc.join
    }


  def likelihoods(root:Root=tree.root):LinearSeq[Double]={
    finalLikelihood(partialLikelihoods(root),m.pi)
  }
  def likelihoods = likelihoods()

  def logLikelihoodRoot(root:Root=tree.root):Double={
    likelihoods(root).zip(aln.countList).map{t=>math.log(t._1)*t._2}.reduceLeft{_+_}
  }
  lazy val logLikelihood:Double={
    logLikelihoodRoot()
  }

  val realCache  = aln.alphabet.matElements.sortBy{_.id}.map{l:Letter=> 
    (0 until numClasses).foldLeft(List.fill(l.alphabet.length * numClasses)(0.0)){(list,i)=>list.updated(l.id +(i * l.alphabet.matLength),1.0)}
  }.toIndexedSeq
  val emptyCache = List.fill(aln.alphabet.length * numClasses)(1.0)

  def leafPartialLikelihoods(l:Letter)={
  l match {
    case a if (a.isReal) => realCache(l.id)
    case a => emptyCache
  }}
  def allLeafPartialLikelihoods(p:List[Pattern],l:Leaf)={val seqId =aln.seqId(l.name); p.map{p2=>leafPartialLikelihoods(p2(seqId))}}

  /*
  val allLeafPartialLikelihoods=immutableHashMapMemo{t:(List[Pattern],Leaf)=> val (p,l)=t
    println("Recalc")
    p.map{p2=> leafPartialLikelihoods(p2(l.id.get))}
  }*/


  def factory(t:Tree,m:Model,aln:Alignment) = new SimpleLikelihoodCalc(t,m,aln)

  def updated(t:Tree)={
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

  def updated(m:Model)={
    factory(tree,m,aln)
  }
}

trait LikelihoodEngine{
  def combinePartialLikelihoods(intermediates:LinearSeq[LinearSeq[PartialLikelihoods]]):LinearSeq[PartialLikelihoods]
  def partialLikelihoodCalc(end:LinearSeq[PartialLikelihoods],matrix:LinearSeq[LinearSeq[Double]]):LinearSeq[PartialLikelihoods]    
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
  def partialLikelihoodCalc(end:LinearSeq[PartialLikelihoods],matrix:LinearSeq[LinearSeq[Double]])={
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
  
  def partialLikelihoodCalc(end:LinearSeq[PartialLikelihoods],matrix:LinearSeq[LinearSeq[Double]])={
    new EnhancedMatrix(end) multipleDotProduct matrix
  }

  def finalLikelihood(partial:LinearSeq[PartialLikelihoods],pi:IndexedSeq[Double])={
    val piList = pi.toList
    partial.map{vec=>
      piList.dotProduct(vec)
    }
  }
}
