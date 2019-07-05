# Lifetable functions

mxax_to_qx <- function(mx, ax){
	qx    <- mx / (1 + (1 - ax) * mx)
	n     <- length(mx)
	qx[n] <- 1
	return(qx)
}

qx_to_px <- function(qx){
	px <- 1 - qx
	return(px)
}

px_to_lx <- function(px){
	n  <- length(px)
	lx <- c(1, cumprod(px[-n]))
	lx
}

mxax_to_lx <- function(mx, ax){
	qx <- mxax_to_qx(mx, ax)
	px <- qx_to_px(qx)
	lx <- px_to_lx(px)
	lx
}

lxqx_to_dx <- function(lx, qx){
	lx * qx
}


lxaxdx_to_Lx <- function(lx, ax, dx){
	Lx <- lx - (1 - ax) * dx
	return(Lx)
}

qxax_to_Lx <- function(qx, ax){
	lx <- px_to_lx(1 - qx)
	dx <- lxqx_to_dx(lx, qx)
	Lx <- lx - (1 - ax) * dx
	return(Lx)
}

Lx_to_Tx <- function(Lx){
	Tx <- rev(cumsum(rev(Lx)))
	return(Tx)
}

Txlx_to_ex <- function(Tx, lx){
	ex <- Tx / lx
	return(ex)
}


my_lt <- function(mx, ax, x){
	qx <- mxax_to_qx(mx, ax)
	px <- qx_to_px(qx)
	lx <- px_to_lx(px)
	dx <- lxqx_to_dx(lx, qx)
	Lx <- qxax_to_Lx(qx, ax)
	Tx <- Lx_to_Tx(Lx)
	ex <- Txlx_to_ex(Tx, lx)
	
	LT <- data.frame(Age = x,
					 mx = mx,
					 ax = ax,
					 qx = qx,
					 px = px,
					 lx = lx,
					 dx = dx,
					 Lx = Lx,
					 Tx = Tx,
					 ex = ex)
	return(LT)
}

my_lt2 <- function(mx, ax, x){
	require(tidyverse)
	data.frame(Age = x, mx = mx, ax = ax) %>% 
		mutate( qx = mxax_to_qx(mx, ax),
				px = qx_to_px(qx),
				lx = px_to_lx(px),
				dx = lxqx_to_dx(lx, qx),
				Lx = qxax_to_Lx(qx, ax),
				Tx = Lx_to_Tx(Lx),
				ex = Txlx_to_ex(Tx, lx))
}

# This one assumes you pass in (LTi) a data.frame
# containing mx and ax at the least
my_lt3 <- function(LTi){
	require(tidyverse)
	LTi %>% 
		mutate( qx = mxax_to_qx(mx, ax),
				px = qx_to_px(qx),
				lx = px_to_lx(px),
				dx = lxqx_to_dx(lx, qx),
				Lx = qxax_to_Lx(qx, ax),
				Tx = Lx_to_Tx(Lx),
				ex = Txlx_to_ex(Tx, lx))
}

# OKOK this isn't a lifetable function,
# but this was easiest
# 
# Lx is literally the Lx column from a lifetable,
# a single vector, BUT the radix of the lifetable
# should be 1!!!
Coales_r <- function(Lx,fxf,x){
	# 1) get R0
	R0 <- sum(Lx * fxf)
	# 2) assume T (but I'll call it G, sorry)
	G   <- 29
	# 3) follow formula for first guess at r:
	ri  <- log(R0) / G
	# 4) now start the loop that updates r:
	for (i in 1:100){
		# step 1:
		delta <- sum(exp(-ri * x) * fxf * Lx) - 1
		# step 2:
		ri <- ri + (delta / (G - (delta / ri)))
	}
	return(ri)
}