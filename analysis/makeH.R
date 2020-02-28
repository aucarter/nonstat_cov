### Define the anistotropy function
makeH <- function(gamma, v_x, v_y) {
    H <- gamma*diag(2) + c(v_x, v_y) %*% t(c(v_x, v_y))
    return(H)
}
# makeH(exp(-0.272), 0.477, -0.313) / exp(-1.75)
