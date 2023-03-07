
#create functions
quart <- function(x) {
  x <- sort(x)
  n <- length(x)
  m <- (n+1)/2
  if (floor(m) != m) {
    l <- m-1/2; u <- m+1/2
  } else {
    l <- m-1; u <- m+1
  }
  c(Q1 = median(x[1:l]), Q2 = median(x), Q3 = median(x[u:n]))
}

normalize.vector <- function(vector) {
  
  if ( near(max(vector), 0) & near(min(vector), 0) ) {
    return(vector)
  } else { 
    return((vector - min(vector))/(max(vector) - min(vector)))
  }
}

my_homogeneity_test <- function(mat) {
  
  groups <- rownames(mat)
  
  for (i in 1:(nrow(mat) - 1)) {
    for (j in (i+1):nrow(mat)) {
      
      dummy.x <- mat[i, ]
      dummy.y <- mat[j, ]
      
      #### sometimes we might not have all bins, but the test is still meaningful
      #### for that particular pair
      nonzero.x <- which(dummy.x > 0)
      nonzero.y <- which(dummy.y > 0)
      
      gd.idcs <- base::union(nonzero.x, nonzero.y)
      dof <- length(gd.idcs) - 1
      
      if (dof > 0) {
        dummy.x <- dummy.x[gd.idcs]
        dummy.y <- dummy.y[gd.idcs]
        
        n1 <- sum(dummy.x)
        n2 <- sum(dummy.y)
        N <- n1 + n2
        
        empirical_p_vec <- (dummy.x + dummy.y) / N
        
        P_1 <- n1 * empirical_p_vec
        Q_1 <- sum( (dummy.x - P_1)^2 / P_1 )
        
        P_2 <- n2 * empirical_p_vec
        Q_2 <- sum( (dummy.y - P_2)^2 / P_2 )
        
        Q <- Q_1 + Q_2
        
        chi.out <- list(p.value = pchisq(Q, df = dof, lower.tail = FALSE))
        
        #### I don't think below is working correctly,
        #### possibly because it is testing for independence instead of homogeneity?
        # chi.out <- chisq.test(x = dummy.x, y = dummy.y, simulate.p.value = TRUE, B = 2000)
        
      } else {
        chi.out <- list(p.value = NA)
      }
      
      
      if (i == 1 & j == 2) {
        dummy_tb <- tibble(group_1 = groups[i], 
                           group_2 = groups[j],
                           p_value = chi.out$p.value)
      } else {
        dummy_tb <- tibble(group_1 = groups[i], 
                           group_2 = groups[j],
                           p_value = chi.out$p.value) %>% 
          bind_rows(dummy_tb)
      }
    }
  }
  
  dummy_tb <- dummy_tb %>% 
    arrange(p_value)
  
  return(dummy_tb)
}

find.pairwise.group.mutation.delta <- function(mat) {
  
  groups <- rownames(mat)
  
  for (i in 1:(nrow(mat) - 1)) {
    for (j in (i+1):nrow(mat)) {
      
      dummy.x <- mat[i, ]
      dummy.y <- mat[j, ]
      
      #### focus on mutations where at least one vector has the mutation
      nonzero.x <- which(dummy.x > 0)
      nonzero.y <- which(dummy.y > 0)
      
      gd.idcs <- base::union(nonzero.x, nonzero.y)
      
      delta.vec <- dummy.x[gd.idcs] - dummy.y[gd.idcs]
      
      if (i == 1 & j == 2) {
        dummy_tb <- tibble(group_1 = groups[i], 
                           group_2 = groups[j],
                           mutation = colnames(mat)[gd.idcs],
                           delta = delta.vec)
      } else {
        dummy_tb <- tibble(group_1 = groups[i], 
                           group_2 = groups[j],
                           mutation = colnames(mat)[gd.idcs],
                           delta = delta.vec) %>% 
          bind_rows(dummy_tb)
      }
    }
  }
  
  dummy_tb <- dummy_tb %>% 
    unite(col = 'comparison', group_1:group_2, sep = '_vs_', remove = FALSE) %>% 
    arrange(abs(delta))
  
  return(dummy_tb)
}
