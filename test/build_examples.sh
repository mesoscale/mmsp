#!/bin/bash

tstart=$(date +%s)
nSerBld=0
nParBld=0
nParRun=0
ITERS=1000
INTER=500
CORES=4
EXEC=true

cd ../examples
examples=$(pwd)

echo -n "Building examples in serial and parallel"


while [[ $# > 0 ]]
do
	key="$1"
	case $key in
		-c|--clean)
		echo -n ", cleaning up after"
		CLEAN=true
		;;
		-p|--purge)
		echo -n ", cleaning up after"
		CLEAN=true
		PURGE=true
		;;
		-n|--noexec)
		echo -n ", not executing"
		EXEC=false
		;;
		-l|--long)
		echo -n ", taking 5,000 steps"
		ITERS=5000
		INTER=1000
		;;
		-x|--extra)
		echo -n ", taking 25,000 steps"
		ITERS=25000
		INTER=5000
		;;
		*)
		# unknown option
		echo "WARNING: Unknown option ${key}."
		echo
    	;;
	esac
	shift # past argument or value
done

echo

exdirs=("coarsening/grain_growth/anisotropic/Monte_Carlo/" \
"coarsening/grain_growth/anisotropic/phase_field/" \
"coarsening/grain_growth/anisotropic/sparsePF/" \
"coarsening/grain_growth/isotropic/Monte_Carlo/" \
"coarsening/grain_growth/isotropic/phase_field/" \
"coarsening/grain_growth/isotropic/sparsePF/" \
"coarsening/ostwald_ripening/isotropic/phase_field/" \
"coarsening/zener_pinning/anisotropic/Monte_Carlo/" \
"coarsening/zener_pinning/anisotropic/phase_field/" \
"coarsening/zener_pinning/anisotropic/sparsePF/" \
"coarsening/zener_pinning/isotropic/Monte_Carlo/" \
"coarsening/zener_pinning/isotropic/phase_field/" \
"coarsening/zener_pinning/isotropic/sparsePF/" \
"phase_transitions/allen-cahn/" \
"phase_transitions/cahn-hilliard/explicit/" \
"phase_transitions/model_A/" \
"phase_transitions/model_B/" \
"phase_transitions/solidification/anisotropic/" \
"phase_transitions/spinodal/" \
"statistical_mechanics/Heisenberg/" \
"statistical_mechanics/Ising/" \
"statistical_mechanics/Potts/")

# Skip phase_transitions/cahn-hilliard/convex_splitting  --  too long
# Skip differential_equations/elliptic/Poisson  --  segmentation faults
# Skip beginners_diffusion  --  does not match execution pattern

n=${#exdirs[@]}
for (( i=0; i<$n; i++ ))
do
	exstart=$(date +%s)
	j=$(($i+1))
	cd $examples/${exdirs[$i]}
	printf "%2d/%2d %-50s\t" $j $n ${exdirs[$i]}

	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	if [[ $EXEC ]]
	then
		rm test.*.png
		mpirun -np $CORES ./parallel --example 2 test.0000.dat && mpirun -np $CORES ./parallel test.0000.dat $ITERS $INTER >/dev/null && ((nParRun++))
		for f in *.dat
		do
			mmsp2png --zoom $f >/dev/null
		done
	fi
	if [[ $CLEAN ]]
	then
		make -s clean
		rm test.*.dat
		if [[ $PURGE ]]
		then
			rm test.*.png
		fi
	fi

	exfin=$(date +%s)
	exlapse=$(echo "$exfin-$exstart" | bc -l)
	echo "${exlapse} seconds"
done
cd ${examples}

echo
echo "${nSerBld} serial   examples built    successfully."
echo "${nParBld} parallel examples built    successfully."
echo "${nParRun} parallel examples executed successfully."

tfinish=$(date +%s)
elapsed=$(echo "$tfinish-$tstart" | bc -l)

echo "${elapsed} seconds elapsed."

cd ../test/
