#!/bin/bash

tstart=$(date +%s)
nSerBld=0
nParBld=0
nParRun=0

cd ../examples

while [[ $# > 0 ]]
do
	key="$1"
	case $key in
		-c|--clean)
		CLEAN=true
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

if [[ $CLEAN ]]; then
	echo "Building examples in serial and parallel, then deleting binaries."
else
	echo "Building examples in serial and parallel, leaving binaries behind (specify --clean to clean up)."
fi
echo

examples=$(pwd)

cd beginners_diffusion/
pwd
	(make -Bs || exit $?) && ((nSerBld++))
	if $CLEAN; then make -s clean; fi
	cd $examples

echo
cd coarsening
pwd
cd grain_growth
pwd
cd anisotropic
pwd
	cd Monte_Carlo
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd ../phase_field/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd ../sparsePF/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

cd coarsening/grain_growth/isotropic
pwd
	cd Monte_Carlo/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd ../phase_field/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd ../sparsePF/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

echo
cd coarsening/ostwald_ripening
pwd
cd isotropic
pwd
	cd phase_field/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

echo
cd coarsening/zener_pinning
pwd
cd anisotropic
pwd
	cd Monte_Carlo/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd ../phase_field/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd ../sparsePF/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

cd coarsening/zener_pinning/isotropic
pwd
	cd Monte_Carlo/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd ../phase_field/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd ../sparsePF/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

echo
cd differential_equations
pwd
cd elliptic
	cd Poisson/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

echo
cd phase_transitions
pwd
	cd allen-cahn/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

	cd phase_transitions/cahn-hilliard
	pwd
	echo "Skipped phase_transitions/cahn-hilliard/convex_splitting/"
# convex splitting is not a suitable example -- it runs much, much too long!
  cd convex_splitting/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
#	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
#	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
#	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd ../explicit/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

	cd phase_transitions/model_A
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

	cd phase_transitions/model_B/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

	cd phase_transitions/spinodal/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

cd phase_transitions/solidification
pwd
	cd anisotropic/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

echo
cd statistical_mechanics
pwd
	cd Heisenberg/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

	cd statistical_mechanics/Ising/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

	cd statistical_mechanics/Potts/
	pwd
	(make -Bs || exit $?) && ((nSerBld++))
	(make -Bs parallel || exit $?) && ((nParBld++))
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && ((nParRun++))
	for f in *.dat; do mmsp2png --zoom $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

echo
echo "${nSerBld} serial   examples built    successfully."
echo "${nParBld} parallel examples built    successfully."
echo "${nParRun} parallel examples executed successfully."

tfinish=$(date +%s)
elapsed=$(echo "$tfinish-$tstart" | bc -l)

echo "${elapsed} seconds elapsed."

cd ../test/
