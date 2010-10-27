import sbt._
import java.io.File



class NewTreeProject(info: ProjectInfo) extends DefaultProject(info) { //with AutoCompilerPlugins {
  override def compileOptions = super.compileOptions ++ Seq(Unchecked,Optimize) 

//  val nativelibs4javaRepo = "NativeLibs4Java Repository" at "http://nativelibs4java.sourceforge.net/maven/"
//  val scalacl = compilerPlugin("com.nativelibs4java" % "scalacl-compiler-plugin" % "1.0-SNAPSHOT") // or "1.0-SNAPSHOT" to get the latest development version

  val scalaToolsSnapshots = "Scala-Tools Maven2 Snapshots Repository" at "http://scala-tools.org/repo-snapshots"
  val clojars = "Clojars Maven2 Repository" at "http://clojars.org/repo/"

  val scalatest = "org.scalatest" % "scalatest" % "1.2-for-scala-2.8.0.RC5-SNAPSHOT"
  val logspace = "Logspace maven repo" at "http://www.logspace.co.uk/maven/"
  val tlf = "uk.co.logspace.tlf" %% "tlf" % "1.2.0"
  val colt = "incanter" % "parallelcolt" % "0.9.4"
  val commonsMath = "org.apache.commons" % "commons-math" % "2.0"
  val commonsCollections = "commons-collections" % "commons-collections" % "3.2.1"
  val beastBeauti = "dr.math" % "beauti" % "1.5.4" from "http://www.logspace.co.uk/jar/beauti-1.5.4.jar"
  val beast = "dr.math" % "beast" % "1.5.4" from "http://www.logspace.co.uk/jar/beast-1.5.4.jar"
  val scalaz = "com.googlecode.scalaz" %% "scalaz-core" % "5.0"

  val concurrent = "concurrent" % "concurrent" % "1.3.4" // needed for proguard
  val jsr166y = "jsr166y" % "jsr166y" % "0" from "http://gee.cs.oswego.edu/dl/jsr166/dist/jsr166y.jar"

  val hawtDispatch = "org.fusesource.hawtdispatch" % "hawtdispatch" % "1.0"

  val jblas = "jblas" % "jblas" % "1.1.1" from "http://github.com/downloads/mikiobraun/jblas/jblas-1.1.1.jar"
 // val sbt = "sbt" % "sbt_2.7.7" % "0.5.7" from "http://simple-build-tool.googlecode.com/svn-history/r1125/artifacts/0.5.7-p1/jars/sbt_2.7.7.jar"
  
}
