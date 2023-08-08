  quad <- function(a, b, c)
  {
    a <- as.complex(a)
    answer <- c((-b + sqrt(b^2 - 4 * a * c)) / (2 * a),
                (-b - sqrt(b^2 - 4 * a * c)) / (2 * a))
    if(all(Im(answer) == 0)) answer <- Re(answer)
    if(answer[1] == answer[2]) return(answer[1])
    if(answer[1] >= 0) return(answer[1])
    if(answer[2]>=0) return(answer[2])
  }
