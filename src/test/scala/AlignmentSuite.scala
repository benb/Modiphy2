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
    val aln = new Fasta(alnStr.lines).toAlignment(AminoAcid)
    aln("Synechocysti").head should be (AminoAcid.G)
    aln.columns.drop(3).head.apply("Gloeobacter") should be (AminoAcid.T)
  }
}
