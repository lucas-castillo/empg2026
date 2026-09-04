compute_gzip <- function(v){
  v <- v[!is.na(v)]
  v <- v - min(v) + 1
  codes <- c((10*16+2):(16*25))
  ascii <- sapply(codes, intToUtf8)
  length(memCompress(paste0(ascii[v], collapse = ""), type = "gzip"))
}

get_Summaries_full <- function(v, ub=NULL, lb=NULL){
  v <- round(v)
  r <- c(FALSE, v[2:length(v)] == v[1:(length(v)-1)])
  if (sum(!r, na.rm=T) == 1){
    return(list(
      "R" = 1,
      "A" = NA,
      "TPF" = NA,
      "D" = NA,
      "S" = NA
    ))
  }
  
  r[1] <- NA
  n <- length(v)
  D <- abs(v[2:n] - v[1:(n-1)])
  A <- (D < (1 + 1e-8)) & (D > (1 - 1e-8))
  A <- c(NA, A)
  D <- c(NA, D)
  A[r] <- NA
  D[r] <- NA
  
  S <- (
    dnorm(v, mean(v, na.rm=T), sd(v, na.rm=T), log=T) -
      dunif(v, min(v, na.rm=T), max(v, na.rm=T), log=T)
  ) / length(v[!is.na(v)])
  
  one <- v[1:(n-2)]
  two <- v[2:(n-1)]
  three <- v[3:(n)]
  
  TP <- (one > two & two < three) | (one < two & two > three)
  TP <- c(NA, NA, TP)
  
  return(list(
    "R" = r,
    "A" = A,
    "TPF" = TP,
    "D" = D,
    "S" = S
  ))
}

.manage_alternatives <- function(a, s){
  # if (!is.null(a) & length(a) == 1){
  #   stop("a should be a vector of alternatives, not the number of alternatives")
  # }
  if (is.null(a)){
    a <- min(s, na.rm=T):max(s, na.rm=T)
  }
  return(a)
}

.response_matrix <- function(s, a, wrap=T){
  M <- matrix(0, length(a), length(a))
  for (i in 2:length(s)){
    last_index <- which(a == s[i-1])
    this_index <- which(a == s[i])
    M[last_index, this_index] <- M[last_index, this_index] + 1
  }
  if (wrap){
    last_index <- which(a == s[i])
    this_index <- which(a == s[1])
    M[last_index, this_index] <- M[last_index, this_index] + 1
  }
  return(M)
}

RNG <- function(s, a=NULL){
  s <- s[!is.na(s)]
  # Evans' RNG. Note T+Neil has a typo here
  a <- .manage_alternatives(a,s)
  M <- .response_matrix(s,a, wrap = T)
  
  row_sums <- rowSums(M)
  top <- sum(M * log10(M), na.rm = T)
  if (top == 0) return(0)
  bottom <- sum(row_sums * log10(row_sums), na.rm = T)
  top / bottom
}

all_measures <- function(s, a=NULL){
  data.frame(t(unlist(
    lapply(
      all_measures_full(s, a), 
      FUN = \(x){apply(x, 2, mean, na.rm=T)}))))
}

all_measures_full <- function(s, a=NULL){
  a <- .manage_alternatives(a,s)
  if (length(a) == 1){
    return(
      list(  
        data.frame(
          R = NA,
          A = NA,
          TPF = NA,
          D = NA,
          S = NA
        ),
        data.frame(
          RNG = NA,
          gzip=NA
        )
      )
    )
  }
  
  summ <- get_Summaries_full(s)
  n <- length(s)
  
  list(  
    data.frame(
      R = summ[[1]],
      A = summ[[2]],
      TPF = summ[[3]],
      D = summ[[4]],
      S = summ[[5]]
    ),
    data.frame(
      RNG = RNG(s, a),
      gzip=compute_gzip(s) / n
    )
  )
}

random_expectations <- function(s, a=NULL){
  ret <- data.frame()
  for (i in 1:1000){
    ret <- rbind(ret, all_measures(sample(s), a))
  }
  data.frame(t(apply(ret, 2, mean)))
}
