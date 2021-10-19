using GaussianRandomFields
using DelimitedFiles
using LinearAlgebra
using DelimitedFiles
using Plots
using NPZ


#-------------------------------------------------
# 	PLOT A SLICE OF THE SAVED ARRAY

function dir__plotslice(dir)
	dims 				=	npzread("$(dir)/dimensions.npy")
	dims 				=	round.(Int64, dims)
	births 				=	npzread("$(dir)/pixel_births.npy")

	println(length(dims))

	if length(dims) == 2
		heatmap(births)
	elseif length(dims) ==3
		heatmap(births[:,:,1])
	end

end

#-------------------------------------------------
# 	FINAL RUN

function dir_pixperaxis_numaxes_covariancetype__save(;
														dir = "/Users/gh10/a/r/ac/re/pu/ar/202006_ffwg_fastforwardalgorithmforgenerators/computation/data/for_paper/cube",
														pixperaxis = 10, 
														numaxes = 2, 
														covariancetype = "exponential"
													)


	# 	GENERATE DATA
	# 	-------------

	dims 					=	fill(pixperaxis, numaxes)

	# 	DETERMINE COVARIANCE FUNCTION
	if covariancetype == "exponential"
		cov 				=	CovarianceFunction(numaxes, Exponential(.5))
	elseif covariancetype == "anisotropic"
		A 					=	fill( 800, (numaxes, numaxes) ) + 200*I 
		cov = CovarianceFunction(numaxes, AnisotropicExponential(A))
	end

	# 	SET PIXEL RANGES
	pts 					=	range(0, stop=1, length=pixperaxis)

	# 	ASSEMBLE GENERATOR
	if numaxes == 2
		grf 				= 	GaussianRandomField(cov, CirculantEmbedding(), pts, pts, minpadding=numaxes*pixperaxis)
	elseif numaxes == 3
		grf 				= 	GaussianRandomField(cov, CirculantEmbedding(), pts, pts, pts, minpadding=numaxes*pixperaxis)
	end
	
	# 	SAMPLE
	field 					=	sample(grf)


	# 	SAVE DATA
	# 	---------	

	# 	CREATE FOLDER
	foldername 				=	"grf_$(covariancetype)$(numaxes)d_$(pixperaxis)pixperaxis"
	dir_save				=	"$(dir)/$(foldername)"
	mkdir(dir_save)

	# 	WRITE CSV
	writedlm("$(dir_save)/dimensions.csv", dims)
	# writedlm("$(dir_save)/pixel_births.csv", arr__rowflattened(field))  # REMEMBER TO FLATTEN DATA

	# 	WRITE NUMPY
	npzwrite("$(dir_save)/pixel_births.npy", field)	
	npzwrite("$(dir_save)/dimensions.npy", dims)		

end



#-------------------------------------------------
# 	FLATTEN ARRAY IN ROW-MAJOR ORDER

function arr__rowflattened(arr)
	dims 					=	size(arr)
	flattened 				=	zeros(prod(dims))

	counter 				=	0

	if length(dims) == 2
		for i0 = 1:dims[1]
			for i1 = 1:dims[2]
				counter 	=	counter+1
				flattened[counter] 	=	arr[i0, i1]
			end
		end

	elseif length(dims) == 3
		for i0 = 1:dims[1]
			for i1 = 1:dims[2]
				for i2 = 1:dims[3]
					counter 	=	counter+1
					flattened[counter] 	=	arr[i0, i1, i2]
				end
			end
		end
	end

	return flattened
end

#-------------------------------------------------
# 	UNFLATTEN ARRAY IN ROW-MAJOR ORDER

function arr_dims__rowUNflattened(arr, dims)
	unflattened				=	zeros(dims...)

	counter 				=	0

	if length(dims) == 2
		for i0 = 1:dims[1]
			for i1 = 1:dims[2]
				counter 	=	counter+1
				unflattened[i0, i1] 	=	arr[counter]
			end
		end

	elseif length(dims) == 3
		for i0 = 1:dims[1]
			for i1 = 1:dims[2]
				for i2 = 1:dims[3]
					counter 	=	counter+1
					unflattened[i0, i1, i2]		=	arr[counter]
				end
			end
		end
	end

	return unflattened
end


