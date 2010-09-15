import org.scalatest.FunSuite
import org.scalatest.matchers.ShouldMatchers._
import modiphy.math._

class MathSuite extends FunSuite {
  test("Identity Matrix"){
    val eye20 = EnhancedIndexedMatrix.eye(20)
    eye20.length should be (20)
    eye20.head.length should be (20)
    eye20(3)(3) should be (1.0)
    eye20(19)(14) should be (0.0)
  }
}
