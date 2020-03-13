### Define the anistotropy function
makeH <- function(gamma, v_x, v_y) {
    H <- gamma*diag(2) + c(v_x, v_y) %*% t(c(v_x, v_y))
    return(H)
}
H <- makeH(exp(-0.272), 0.477, -5) / exp(-1.75)

plotH <- function(H) {
    eigenH <- eigen(H)
    eigVal <- eigenH$values
    eigVec <- eigenH$vectors
    eigScl  <- eigVec %*% diag(sqrt(eigVal))  # scale eigenvectors to length = square-root
    ctr    <- c(0, 0)
    angles <- seq(0, 2*pi, length.out=200)
    xMat    <- rbind(ctr[1] + eigScl[1, ], ctr[1] - eigScl[1, ])
    yMat    <- rbind(ctr[2] + eigScl[2, ], ctr[2] - eigScl[2, ])
    ellBase <- cbind(sqrt(eigVal[1])*cos(angles), sqrt(eigVal[2])*sin(angles)) # normal ellipse
    ellRot  <- eigVec %*% t(ellBase)                                          # rotated ellipse
    plot((ellRot+ctr)[1, ], (ellRot+ctr)[2, ], asp=1, type="l", lwd=2, xlab = "x", ylab = "y")
    matlines(xMat, yMat, lty=1, lwd=2, col="green")
    points(ctr[1], ctr[2], pch=4, col="red", lwd=3)
}

