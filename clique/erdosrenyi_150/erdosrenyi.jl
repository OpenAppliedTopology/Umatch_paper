
using Random

function generate_erdosrenyi(n)

	dismat 			=	zeros(n,n)
	permlength 		=	(n^2 - n)/2
	permlength 		=	round(Int64, permlength)
	perm 			=	randperm(permlength)

	counter 		=	0
	for i = 1:n
		for j = (i+1):n
			counter =	counter+1
			dismat[i,j] 	=	perm[counter]
			dismat[j,i] 	=	perm[counter]
		end
	end

	return dismat
	
end