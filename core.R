
# Maps values to interpolated color gradient
colorGradient = function(x,
	colors=c("red","yellow","green"),
	colsteps=100)
{
	return(
		colorRampPalette(colors)(colsteps) [
			findInterval(x, seq(min(x), max(x), length.out=colsteps))
		]
	)
}

# Mix vectors of two colors by minimum threshold method
mixColors = function(colors1, colors2) {
	col_mat1 = col2rgb(colors1)
	col_mat2 = col2rgb(colors2)

	col_mat_combined = pmin(col_mat1, col_mat2)

	return(apply(col_mat_combined, 2, function(x) rgb(x[1], x[2], x[3], maxColorValue=255)))
}

# Clamp negative values to minimum observed non-zero value per variable
clampZero = function(mat) {

	# Get lowest non-zero value per column
	min_non_neg = apply(mat, 2, function(col) min(col[col > 0]))

	# Clamp matrix by column
	for (i in 1:ncol(mat)) {
		zero_idx = mat[, i] == 0  # finds netagive values
		mat[zero_idx, i] = min_non_neg[i]  # replace with non-zero value
	}

	return(mat)
}


# Scatterplot of double positives, interpretted as errors
# Estimates error rate
plotDoublePositive = function(x, y,
	xcut, ycut,
	error_color=brewer.pal(9, "Set1")[1],
	...
) {
	if (length(x) != length(y)) {
		stop("x and y vectors have different dimensions.")
	}
	
	error = x > xcut & y > ycut
	pt_col = rep("black", length(x))
	pt_col[error] = error_color

	pt_pch = rep(1, length(x))
	pt_pch[error] = 16

	plot(x, y,
		col=pt_col,
		pch=pt_pch,
		...
	)

	grey = rgb(100, 100, 100, maxColorValue=255)
	abline(v=xcut, lty=2, col=grey)
	abline(h=ycut, lty=2, col=grey)

	error_rate = mean(error)

	legend("topright",
		legend=paste(format(error_rate * 100, digits=3), "%"),
		text.col=error_color,
		bty="n")
}


# Plots 2 pairwise scatter plots
# First column is assumed to be H2B stain
pairwiseNullPlot = function(mat,
	cutoffs, ...
) {
	if (ncol(mat) != 4) {
		stop("Matrix dimension is not 4")
	}

	par(mfrow=c(1, 3))

	for (coord in list(
		c(3, 2), c(4, 2), c(3, 4)))
	{
		i = coord[1]
		j = coord[2]
		plotDoublePositive(x=mat[, i], y=mat[, j],
			xcut=cutoffs[i],
			ycut=cutoffs[j],
			xlab=colnames(mat)[i],
			ylab=colnames(mat)[j],
			...
		)
	}
}

