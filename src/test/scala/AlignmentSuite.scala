import org.scalatest.FunSuite
import org.scalatest.matchers.ShouldMatchers._
import modiphy.alignment._

class AlignmentSuite extends FunSuite {

  test("DNA"){
    DNA.length should be (4)
    DNA.A should not be (DNA.C)
    DNA.A should be (DNA.A)
  }
  test("Protein"){
    AminoAcid.length should be (20)
    AminoAcid.C should not be (DNA.C)
  }
  test("Fasta"){
    val alnStr =  """>Chlamydomona
    GRIVQIIGPVLDIVF
    >Chlorella
    GRITQIIGPVLDVSF
    >Cyanidioschy
    GKVVQIIGPVLDIEF
    >Cyanidium     
    GKVIQIIGPVLDIVF
    >Ecoli        
    GKIVQVIGAVVDVEF
    >Euglena
    GIILQIIGPVMDISF
    >Gloeobacter
    GYITQIIGPVVDAEF
    >Guillardia
    GYITQIIGPVVDVEF
    >Mycobacteriu
    GRVVRVTGPVVDVEF
    >Mycoplasma 
    GKVHQVIGPVVDVIF
    >Nephroselmis 
    GKIVQIVGPVMDVAF
    >Porphyra  
    GSVTQIIGPVLDIAF
    >Prochlorococ
    GVIRQVIGPVLDVEF
    >Rickettsia 
    GKITQVISAVVDVKF
    >Synechocysti 
    GKITQVIGPVIDAQF
    >Thermosynech
    GFITQVIGPVVDIEF
    """
    val aln = new Fasta(alnStr.lines).parseWith(AminoAcid)
    aln("Synechocysti").head should be (AminoAcid.G)
    aln.columns.drop(3).head.apply(aln.seqId("Gloeobacter")) should be (AminoAcid.T)
    aln.frequencies.reduceLeft(_+_) should be (1.0 plusOrMinus 1E-7)
    aln.restrictTo(List("Cyanidium","Ecoli","Gloeobacter")).numSeqs should be(3)
  }
}
