#!/bin/bash

# build_examples.sh
# This test script compiles each canonical* example in serial and parallel
# then executes the example program in 2D. The length of the test can be
# controlled by selecting default (1,000 steps), long (5,000 steps), or extra
# long (25,000 steps). The resulting checkpoints are then converted to PNG
# images to spot check problems in the code or environment, including execution
# errors and numerical instabilities. If you see something, say something.
#
# Questions/comments to trevor.keller@gmail.com (Trevor Keller)

# * Canonical examples exclude the following directories:
# differential_equations/elliptic/Poisson  --  segmentation faults
# beginners_diffusion  --  does not match execution pattern

# Valid flags are:
# --noexec  build in serial and parallel, but do not execute
# --noviz   do not convert data for visualization
# --force   pass -B flag to make
# --np X    pass mpirun X ranks (e.g., --np 3 yields mpirun -np 3)
# --short   execute tests .1x default
# --long    execute tests  5x longer  than default
# --extra   execute tests 25x longer  than default
# --clean   delete binary files after test completes
# --purge   delete binaries, data, and images after test completes

# Initialize timer and completion counters
tstart=$(date +%s)
nSerErr=0
nParErr=0
nRunErr=0
nSerBld=0
nParBld=0
nParRun=0
MFLAG="-s"

# Set execution parameters
ITERS=1000
INTER=500
CORES=4
COREMAX=$(nproc)
if [[ $CORES > $COREMAX ]]
then
	CORES=$COREMAX
fi

# Get going
cd ../examples
examples=$(pwd)

echo -n "Building examples in serial and parallel"


while [[ $# > 0 ]]
do
	key="$1"
	case $key in
		--force)
			echo -n ", forcing build"
			MFLAG="-Bs"
		;;
			--clean)
			echo -n ", cleaning up after"
		CLEAN=true
		;;
		--purge)
			echo -n ", cleaning up after"
			CLEAN=true
			PURGE=true
		;;
		--noexec)
			echo -n ", not executing"
			NEXEC=true
		;;
		--short)
			ITERS=$(($ITERS/10))
			INTER=100
		;;
		--long)
			ITERS=$((5*$ITERS))
			INTER=1000
		;;
		--extra)
			ITERS=$((25*$ITERS))
			INTER=5000
		;;
		--noviz)
			echo -n ", no PNG output"
			NOVIZ=true
		;;
		--np)
			shift
			CORES=$1
		;;
		*)
			echo "WARNING: Unknown option ${key}."
			echo
    ;;
	esac
	shift # pop first entry from command-line argument list, reduce $# by 1
done

if [[ ! $NEXEC ]]
then
		echo ", taking $ITERS steps, using $CORES/$COREMAX MPI ranks"
else
	echo
fi

if [[ ! $NOVIZ ]]
then
	# Remove existing images first. If the test fails,
	# no images in the directory is an obvious red flag.
	rm -f test.*.png
	if [[ $(which mmsp2png) == "" ]]
	then
		# Consult doc/MMSP.manual.pdf if this fails.
		echo "mmsp2png utility not found. Please check your installation,"
		echo " or pass --noviz flag to suppress PNG output."
	fi
fi

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
"phase_transitions/cahn-hilliard/convex_splitting/" \
"phase_transitions/cahn-hilliard/explicit/" \
"phase_transitions/model_A/" \
"phase_transitions/model_B/" \
"phase_transitions/solidification/anisotropic/" \
"phase_transitions/spinodal/" \
"statistical_mechanics/Heisenberg/" \
"statistical_mechanics/Ising/" \
"statistical_mechanics/Potts/")

n=${#exdirs[@]}
for (( i=0; i<$n; i++ ))
do
	exstart=$(date +%s)
	j=$(($i+1))
	cd $examples/${exdirs[$i]}
	printf "%2d/%2d %-52s\t" $j $n ${exdirs[$i]}
	echo $(date) >test.log
	if make $MFLAG
	then
		((nSerBld++))
	else
		((nSerErr++))
	fi
	if make $MFLAG parallel
	then
		((nParBld++))
	else
		((nParErr++))
	fi
	if [[ -f parallel ]] && [[ ! $NEXEC ]]
	then
		# Run the example in parallel, for speed.
		mpirun -np $CORES ./parallel --example 2 test.0000.dat 1>test.log 2>error.log \
		&& mpirun -np $CORES ./parallel test.0000.dat $ITERS $INTER 1>>test.log 2>>error.log
		# Return codes are not reliable. Litter the hard drive with files instead.
		if [[ -f error.log ]] && [[ $(wc -w error.log) > 1 ]]
		then
			((nRunErr++))
			wc -w error.log
			if [[ -f error.log ]]
			then
				head error.log
			fi
		else
			((nParRun++))
			rm -f error.log
		fi
		if [[ ! $NOVIZ ]]
		then
			# Show the result
			for f in *.dat
			do
				mmsp2png --zoom $f >>test.log
			done
		fi
	fi
	# Clean up binaries and images
	if [[ $CLEAN ]]
	then
		make -s clean
		rm -f test.*.dat
		rm -f test.log
		rm -f error.log
		if [[ $PURGE ]]
		then
			rm -f test.*.png
		fi
	fi

	exfin=$(date +%s)
	exlapse=$(echo "$exfin-$exstart" | bc -l)
	printf "%3d seconds\n" $exlapse
done

cd ${examples}

tfinish=$(date +%s)
elapsed=$(echo "$tfinish-$tstart" | bc -l)

printf "Elapsed time: %53d seconds\n" $elapsed
echo
printf "%2d serial   examples compiled successfully, %2d failed.\n" $nSerBld $nSerErr
printf "%2d parallel examples compiled successfully, %2d failed.\n" $nParBld $nParErr
printf "%2d parallel examples executed successfully, %2d failed.\n" $nParRun $nRunErr
cd ../test/

AllERR=$(echo "$nSerErr+$nParErr+$nRunErr" | bc -l)
if [[ $AllERR > 0 ]]
then
	echo "${AllERR} tests failed."
	echo
	exit 1
fi
