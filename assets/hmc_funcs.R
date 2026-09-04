gradient <- function(f, x){pracma::grad(f, x)}
dotProduct <- function(a,b){a %*% b}
joint_density <- function(theta, momentum, log_func){
  log_func(theta) - (.5 * dotProduct(momentum, momentum))
}

leapfrog_step <- function(theta, momentum, epsilon, L, log_pdf){
  thetas <- matrix(nrow=L+1, ncol=length(theta))
  thetas[1,] <- theta
  # start with half step for momentum
  momentum = momentum + (epsilon/2) * gradient(log_pdf, theta);
  # alternate full steps for position and momentum
  for (i in 1:L){
    theta = theta + epsilon * momentum
    thetas[i+1,] <- theta
    if (i != L){
      # full step for momentum except for the end of trajectory
      momentum = momentum + epsilon * gradient(log_pdf, theta)
    }
  }
  # make a half step (instead of a full one) for the momentum at the end
  momentum = momentum + (epsilon/2) * gradient(log_pdf, theta)
  
  # negate momentum to make the proposal symmetric
  momentum = -1 * momentum;
  return(thetas)
}
