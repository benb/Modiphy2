package modiphy.tree
import scala.collection.immutable.Map.WithDefault
import modiphy.alignment._
import modiphy.model._
import modiphy.alignment.GlobalAlphabet._
import scala.collection.LinearSeq

import modiphy.opt._

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

class MixtureLikelihoodCalc(tree:Tree,aln:Alignment,m:Model,lkl:Option[Seq[LikelihoodCalc[Model]]]=None) extends LikelihoodCalc[Model]{
  def priors = m.priors
  val models = m.models
  val lklCalc = lkl.getOrElse{models.map{new SimpleLikelihoodCalc(tree,_,aln)}}
  
  lazy val likelihoods  = {
    val myLikelihoods = if (Parallel.on && false){
      import jsr166y._
      class Calc(subModels:List[(LikelihoodCalc[Model],Double)]) extends RecursiveTask[Seq[Iterator[Double]]]{
        def compute = {
          subModels match {
            case model::Nil => subModels.map{t=> t._1.likelihoods.map{_ * t._2}}.map{_.iterator}
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
      }.map{_.iterator}
    }
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

/*
trait CachingLikelihoodCalc[A <: Model] extends LikelihoodCalc[A]{
   val partialLikelihoodCache:modiphy.util.Memo[RootedTreePosition,LinearSeq[PartialLikelihoods]]=new modiphy.util.ArrayMemo[RootedTreePosition,LinearSeq[PartialLikelihoods]]({treePos=>treePos.id},tree.maxID)({ treePos => realPartialLikelihoodCalc(treePos)})
   override def partialLikelihoods(treePos:RootedTreePosition)=partialLikelihoodCache(treePos)  
   def tree:Tree
   def realPartialLikelihoodCalc(treePos:RootedTreePosition):LinearSeq[PartialLikelihoods]
}*/
trait LikelihoodCalc[A <: Model]{
//   def partialLikelihoods(treePos:RootedTreePosition):LinearSeq[PartialLikelihoods]
   def updated(t:Tree):LikelihoodCalc[A]
   def updated(m:A):LikelihoodCalc[A]
   def logLikelihood:Double
   def updatedVec(p:VectorParamName,vec:IndexedSeq[Double],paramIndex:ParamMatcher):LikelihoodCalc[A]
   def updatedMat(p:MatrixParamName,mat:IndexedSeq[IndexedSeq[Double]],paramIndex:ParamMatcher):LikelihoodCalc[A]
   def updatedSingle(p:SingleParamName,d:Double,paramIndex:ParamMatcher):LikelihoodCalc[A]
   def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:ParamMatcher):LikelihoodCalc[A]
   def model:A
   def likelihoods:LinearSeq[Double]
}
class SimpleLikelihoodCalc(val tree:Tree,m:Model,val aln:Alignment,val engine:LikelihoodEngine=DefaultLikelihoodFactory.apply) extends LikelihoodCalc[Model]{
  import SimpleLikelihoodCalc._
  
  import engine.combinePartialLikelihoods
  import engine.partialLikelihoodCalc
  import engine.finalLikelihood

  val numClasses = m.numClasses

  
  def model=m
  def setOptParam(p:ParamName,vec:IndexedSeq[Double],paramIndex:ParamMatcher)={ updated(m.setOptParam(p,vec,paramIndex))}

   def updatedVec(p:VectorParamName,vec:IndexedSeq[Double],paramIndex:ParamMatcher)={ updated(m.updatedVec(p,vec,paramIndex)) }
   def updatedSingle(p:SingleParamName,d:Double,paramIndex:ParamMatcher)={ updated(m.updatedSingle(p,d,paramIndex)) }
   def updatedMat(p:MatrixParamName,mat:IndexedSeq[IndexedSeq[Double]],paramIndex:ParamMatcher)={ updated(m.updatedMat(p,mat,paramIndex)) }
  def realPartialLikelihoodCalc(treePos:RootedTreePosition):LinearSeq[PartialLikelihoods]={
    val p = aln.patternList
       (treePos,treePos.get) match {
          case (myTP:TreePositionDir,n:INode)=>
           partialLikelihoodCalc(
            combinePartialLikelihoods(treePos.children.toList.map{realPartialLikelihoodCalc}
           ),m(myTP.upEdge))
          case (myTP:RootedTreePosition,n:INode)=>
            combinePartialLikelihoods(treePos.children.toList.map{realPartialLikelihoodCalc})
          case (myTP:TreePositionDir,l:Leaf)=>   
          partialLikelihoodCalc(
            p.map{p2=>
              leafPartialLikelihoods(p2(l.id.get))},
              m(myTP.upEdge)
            )
        }
    }



    import scalaz._
    import Scalaz._

//   val partialLikelihoods = mutableHashMapMemo(realPartialLikelihoodCalc)
//   val partialLikelihoods:modiphy.util.Memo[RootedTreePosition,LinearSeq[PartialLikelihoods]]=modiphy.util.Memo[RootedTreePosition,LinearSeq[PartialLikelihoods]]({ treePos => realPartialLikelihoodCalc(treePos)})
//   val partialLikelihoods:modiphy.util.Memo[RootedTreePosition,LinearSeq[PartialLikelihoods]]=new modiphy.util.ArrayMemo[RootedTreePosition,LinearSeq[PartialLikelihoods]]({treePos=>treePos.id},tree.maxID)({ treePos => realPartialLikelihoodCalc(treePos)})
    
    def partialLikelihoods(treePos:RootedTreePosition)={
      Parallel.on match {
        case true => parallelPartialLikelihoods(treePos)
        case false => realPartialLikelihoodCalc(treePos)
      }
    }
    def parallelPartialLikelihoods(treePos:RootedTreePosition)={
      import jsr166y._
      import Parallel._
      val p = aln.patternList
      class Calc(treePos:RootedTreePosition) extends RecursiveTask[LinearSeq[PartialLikelihoods]]{
        def compute():LinearSeq[PartialLikelihoods]={
          (treePos,treePos.get) match {
            case (tp:TreePositionDir,n:INode)=>
                val childCalc = treePos.children.toList.map{new Calc(_)}.map{_.fork}
                partialLikelihoodCalc( combinePartialLikelihoods(childCalc.map{_.join}) , m(tp.upEdge) )
            case (tp:RootedTreePosition,n:INode)=>
              val childCalc = treePos.children.toList.map{new Calc(_)}.map{_.fork}
              combinePartialLikelihoods( childCalc.map{_.join})
            case (tp:TreePositionDir,l:Leaf)=> partialLikelihoodCalc(allLeafPartialLikelihoods(p,l),m(tp.upEdge))
          }
        }
      }
      val calc = new Calc(treePos)
      forkJoinPool submit calc
      calc.join
    }


  def likelihoods(root:RootedTreePosition=tree.traverseDown(tree.defaultRoot)):LinearSeq[Double]={
    finalLikelihood(partialLikelihoods(root),m.pi(root.get))
  }
  def likelihoods = likelihoods()

  def logLikelihoodRoot(root:RootedTreePosition=tree.traverseDown(tree.defaultRoot)):Double={
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
  def allLeafPartialLikelihoods(p:List[Pattern],l:Leaf)=p.map{p2=>leafPartialLikelihoods(p2(l.id.get))}

  /*
  val allLeafPartialLikelihoods=immutableHashMapMemo{t:(List[Pattern],Leaf)=> val (p,l)=t
    println("Recalc")
    p.map{p2=> leafPartialLikelihoods(p2(l.id.get))}
  }*/


  def factory(t:modiphy.tree.Tree,m:Model,aln:Alignment) = new SimpleLikelihoodCalc(t,m,aln)

  def updated(t:modiphy.tree.Tree)={
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
    factory(tree,m,aln)//,Map[RootedTreePosition,Map[Seq[Letter],PartialLikelihoods]]())
  }
}

trait LikelihoodEngine{
  def combinePartialLikelihoods(intermediates:LinearSeq[LinearSeq[PartialLikelihoods]]):LinearSeq[PartialLikelihoods]
  def partialLikelihoodCalc(end:LinearSeq[PartialLikelihoods],matrix:Matrix):LinearSeq[PartialLikelihoods]    
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
  def partialLikelihoodCalc(end:LinearSeq[PartialLikelihoods],matrix:Matrix)={
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
  
  def partialLikelihoodCalc(end:LinearSeq[PartialLikelihoods],matrix:Matrix)={
    new EnhancedMatrix(end) multipleDotProduct matrix
  }

  def finalLikelihood(partial:LinearSeq[PartialLikelihoods],pi:IndexedSeq[Double])={
    val piList = pi.toList
    partial.map{vec=>
      piList.dotProduct(vec)
    }
  }
}
