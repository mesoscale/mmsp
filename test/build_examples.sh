#!/bin/bash

tstart=$(date +%s)

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
	make -Bs || exit $?
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
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd ../phase_field/
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd ../sparsePF/
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

cd coarsening/grain_growth/isotropic
pwd
	cd Monte_Carlo/
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd ../phase_field/
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd ../sparsePF/
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
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
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
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
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd ../phase_field/
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd ../sparsePF/
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

cd coarsening/zener_pinning/isotropic
pwd
	cd Monte_Carlo/
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd ../phase_field/
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd ../sparsePF/
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

echo
cd differential_equations
pwd
cd elliptic
pwd
	cd Poisson/
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

echo
cd phase_transitions
pwd
	cd allen-cahn/
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

	cd phase_transitions/cahn-hilliard
	pwd
	cd convex_splitting/
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd ../explicit/
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

	cd phase_transitions/model_A
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

	cd phase_transitions/model_B/
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

	cd phase_transitions/spinodal/
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

cd phase_transitions/solidification
pwd
	cd anisotropic/
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

echo
cd statistical_mechanics
pwd
	cd Heisenberg/
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

	cd statistical_mechanics/Ising/
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

	cd statistical_mechanics/Potts/
	pwd
	make -Bs || exit $?
	make -Bs parallel || exit $?
	mpirun -np 4 ./parallel --example 2 test.0000.dat && mpirun -np 4 ./parallel test.0000.dat 1000 500 >/dev/null && for f in *.dat; do mmsp2png $f >/dev/null; done
	rm test.*.dat
	if $CLEAN; then make -s clean; fi
	cd $examples

tfinish=$(date +%s)
elapsed=$(echo "$tfinish-$tstart" | bc -l)

echo
echo "Examples built successfully. ${elapsed} seconds elapsed."

cd ../test/
