#!/bin/bash
#rm -r __pycache__/
#python traceplt.py --workDir $1 				&
#python statplt.py --workDir $1	 				&
#./dynamicsMovie.sh $1 						&
#./continuityMovie.sh $1 					&
#./innovationMovie.sh $1 					&

~/julia-1.7.1/bin/julia --project=. tracePlot.jl $1 		&
~/julia-1.7.1/bin/julia --project=. errorPlot.jl $1 		&
~/julia-1.7.1/bin/julia --project=. dynamicsMovie.jl $1 	&
~/julia-1.7.1/bin/julia --project=. continuityMovie.jl $1 	&
~/julia-1.7.1/bin/julia --project=. innovationMovie.jl $1 	&
