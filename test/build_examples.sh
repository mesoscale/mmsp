#!/bin/bash

tstart=$(date +%s)
nSerBld=0
nParBld=0
nParRun=0
iters=1000

cd ../examples

while [[ $# > 0 ]]
do
	key="$1"
	case $key in
		-c|--clean)
		CLEAN=true
		shift # past argument
		;;
		-l|--long)
		iters=5000
		shift # past argument
		;;
		*)
		# unknown option
		echo "WARNING: Unknown option ${key}."
		echo
    	;;
	esac
	shift # past argument or value
done

if [[ $CLEAN ]]
then
	echo "Building examples in serial and parallel, then deleting binaries."
else
	echo "Building examples in serial and parallel, leaving binaries behind (specify --clean to clean up)."
fi
echo

examples=$(pwd)

exdirs="coarsening/grain_growth/anisotropic/Monte_Carlo/
coarsening/grain_growth/anisotropic/phase_field/
coarsening/grain_growth/anisotropic/sparsePF/
coarsening/grain_growth/isotropic/Monte_Carlo/
coarsening/grain_growth/isotropic/phase_field/
coarsening/grain_growth/isotropic/sparsePF/
coarsening/ostwald_ripening/isotropic/phase_field/
coarsening/zener_pinning/anisotropic/Monte_Carlo/
coarsening/zener_pinning/anisotropic/phase_field/
coarsening/zener_pinning/anisotropic/sparsePF/
coarsening/zener_pinning/isotropic/Monte_Carlo/
coarsening/zener_pinning/isotropic/phase_field/
coarsening/zener_pinning/isotropic/sparsePF/
phase_transitions/allen-cahn/
phase_transitions/cahn-hilliard/explicit/
phase_transitions/model_A/
phase_transitions/model_B/
phase_transitions/solidification/anisotropic/
phase_transitions/spinodal/
statistical_mechanics/Heisenberg/
statistical_mechanics/Ising/
statistical_mechanics/Potts/"

# Skip phase_transitions/cahn-hilliard/convex_splitting  --  too long
# Skip differential_equations/elliptic/Poisson  --  segmentation faults
# Skip beginners_diffusion  --  does not match execution pattern

for d in $exdirs
do
	cd ${examples}/$d
	printf "%-55s\t" $d
	exstart=$(date +%s)

	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat $iters 500 >/dev/null && ((nParRun++))

	for f in *.dat
	do
		mmsp2png --zoom $f >/dev/null
	done

	if [[ $CLEAN ]]
	then
		rm test.*.dat
		rm test.*.png
		make -s clean
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
