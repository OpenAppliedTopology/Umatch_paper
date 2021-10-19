
using Distances

function get_3dtorus_distance_from_cols(X)
	
	dismat 		=	pairwise(Euclidean(), X)

	for i = [-1, 0, 1]
		for j = [-1, 0, 1]
			for k = [-1, 0, 1]
				Y 				=	X .+ [i,j,k]
				dismat_new		=	pairwise(Euclidean(), Y)
				dismat 			=	min.(dismat, dismat_new)
			end
		end
	end

	return dismat

end