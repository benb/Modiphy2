package modiphy.util

class Memo[A, B](f:A=>B) extends Function[A,B]{
  var map=Map[A,B]()
  def apply(a:A)={
    if (map contains a){map(a)}
    else {val ans = f(a);map=map updated (a,ans);ans}
  }
}
