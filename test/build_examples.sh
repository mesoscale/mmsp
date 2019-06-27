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
# beginners_diffusion  --  does not match execution pattern
# differential_equations/elliptic/Poisson  --  serial execution only
# These will be compiled, but not executed.

# Valid flags are:
# --noexec  build in serial and parallel, but do not execute
# --noviz   do not convert data for visualization
# --pvd     convert to PVD and VTI instead of PNG
# --force   pass -B flag to make
# --np X    pass mpirun X ranks (e.g., --np 3 yields mpirun -np 3)
# --1D      initialize and run examples in 1D
# --2D      initialize and run examples in 2D
# --3D      initialize and run examples in 3D
# --short   execute tests .1x default
# --long    execute tests  5x longer  than default
# --extra   execute tests 25x longer  than default
# --clean   delete generated files (binaries, data, imges) after test completes

# Set output colors
RED='\033[0;31m'  # red
BRED='\033[1;31m'  # bold red
GRN='\033[0;32m'  # green
BGRN='\033[1;32m'  # bold green
YLW='\033[0;33m'  # yellow
BYLW='\033[1;33m' # bold yellow
WHT='\033[0m'     # normal

# Initialize timer and completion counters
tstart=$(date +%s)
nSerErr=0
nParErr=0
nRunErr=0
nSerBld=0
nParBld=0
nParRun=0
MFLAG="-s"

# Default execution parameters
DIM=2            # grid dimensions
ITERS=1000       # number of steps
INTER=500        # steps between checkpoint
ZED="0000"       # checkpoint filename format
CORES=4          # MPI ranks to use
COREMAX=$(nproc) # MPI ranks available
if [[ $CORES -gt $COREMAX ]]
then
	# A lot of machines, Travis-CI build bots included,
	# have only 2 cores available. Let's be polite.
	CORES=$COREMAX
fi

# Get going
root_dir=$(pwd)/../

echo "Building utilities"
cd ${root_dir}/utility
make

cd ${root_dir}/examples
examples=$(pwd)

echo -n "Building examples in serial and parallel"

while [[ $# -gt 0 ]]
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
		--noexec)
			echo -n ", not executing"
			NEXEC=true
		;;
		--1D)
			DIM=1
		;;
		--2D)
			DIM=2
		;;
		--3D)
			DIM=3
		;;
		--np)
			shift
			CORES=$1
		;;
		--short)
			ITERS=$(($ITERS/10))
			INTER=100
			ZED="000"
		;;
		--long)
			ITERS=$((5*$ITERS))
			INTER=1000
		;;
		--extra)
			ITERS=$((25*$ITERS))
			INTER=1000
			ZED="00000"
		;;
		--noviz)
			echo -n ", no PNG output"
			NOVIZ=true
		;;
		--pvd)
			PVD=true
		;;
		*)
			echo -ne ", ${BYLW}ignoring unknown option ${key}${WHT}"
    ;;
	esac
	shift # pop first entry from command-line argument list, reduce $# by 1
done

if [[ ! $NEXEC ]]
then
	echo ", ${DIM}-dimensional, taking $ITERS steps, using $CORES/$COREMAX MPI ranks"
	if [[ $DIM -eq 3 ]] && [[ $ITERS -gt 100 ]]
	then
		DENOM=$((6*$CORES))
		if [[ $ITERS -gt 10000 ]]
		then
			# Not as much work at long times
			DENOM=$((9*$CORES))
		fi
		RTIM=$(($ITERS/$DENOM))
		echo -e "${BYLW}      Expected runtime is $RTIM minutes. Consider --short.${WHT}"
	fi
else
	echo
fi

if [[ ! $NOVIZ ]]
then
	if [[ ! $PVD ]] && [[ ! -f ${root_dir}/utility/mmsp2png ]]
	then
		# Consult doc/MMSP.manual.pdf if this fails.
		echo -e "${BYLW}mmsp2png utility not found. Please check your installation, or pass --noviz.${WHT}"
	fi
	if [[ $PVD ]] && [[ ! -f ${root_dir}/utility/mmsp2pvd ]]
	then
		# Consult doc/MMSP.manual.pdf if this fails.
		echo -e "${BYLW}mmsp2pvd utility not found. Please check your installation, or pass --noviz.${WHT}"
	fi
	if [[ $DIM -eq 3 ]] && [[ ! $PVD ]]
	then
		echo -e "${BYLW}Pass --pvd for ParaView Data output, recommended over PNG for 3D data.${WHT}"
	fi
fi

echo "---------------------------------------------------------------------------"

exdirs=("beginners_diffusion/" \
"HeatConduction/heterogeneous" \
"HeatConduction/homogeneous" \
"coarsening/grain_growth/anisotropic/Monte_Carlo/" \
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
"differential_equations/elliptic/Poisson/" \
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
	if [[ ${exdirs[$i]} != *"beginners"* ]]
	then
		if make $MFLAG parallel
		then
			((nParBld++))
		else
			((nParErr++))
		fi
	fi
	if [[ -f parallel ]] && [[ ! $NEXEC ]]
	then
		rm -f test.*.dat test.*.png test.*.vti test.pvd
		RUNCMD="mpirun -np $CORES ./parallel"
		if [[ ${exdirs[$i]} == *"Poisson"* ]]
		then
			RUNCMD="./poisson"
		fi
		# Run the example in parallel, for speed.
		$RUNCMD --example $DIM test.$ZED.dat 1>test.log 2>error.log \
		&& $RUNCMD test.$ZED.dat $ITERS $INTER 1>>test.log 2>>error.log
		# Return codes are not reliable. Save errors to disk for postmortem.
		if [[ -f error.log ]] && [[ $(wc -w error.log) > 1 ]]
		then
			echo -e "${RED} --FAILED--${WHT}"
			((nRunErr++))
			if [[ -f error.log ]]
			then
				echo "      error.log has the details (head follows)"
				head error.log | sed -e 's/^/      /'
			fi
		else
			((nParRun++))
			rm -f error.log
			if [[ ! $NOVIZ ]]
			then
				# Show the result
				for f in *.dat
				do
					${root_dir}/utility/./mmsp2png --zoom $f >>test.log
				done
				if [[ $PVD ]]
				then
					${root_dir}/utility/./mmsp2pvd --output=test.pvd test.*.dat >>test.log
				fi
			fi
			exfin=$(date +%s)
			exlapse=$(echo "$exfin-$exstart" | bc -l)
			printf "${GRN}%3d seconds${WHT}\n" $exlapse
		fi
	elif [[ ${exdirs[$i]} == *"beginners"* ]]
	then
		exfin=$(date +%s)
		exlapse=$(echo "$exfin-$exstart" | bc -l)
		if [[ ! $NEXEC ]]
		then
			printf "${YLW}%3d seconds\n      example does not match generic pattern for execution${WHT}\n" $exlapse
		else
			printf "${GRN}%3d seconds${WHT}\n" $exlapse
		fi
	else
		exfin=$(date +%s)
		exlapse=$(echo "$exfin-$exstart" | bc -l)
		printf "${GRN}%3d seconds${WHT}\n" $exlapse
	fi
	# Clean up binaries and images
	if [[ $CLEAN ]]
	then
		make -s clean
		rm -f test.*.dat
		rm -f test.log error.log
		rm -f test.*.png test.pvd test.*.vti
	fi
done

cd ${examples}

tfinish=$(date +%s)
elapsed=$(echo "$tfinish-$tstart" | bc -l)
echo "---------------------------------------------------------------------------"
printf "Elapsed time: %53d seconds\n" $elapsed
echo
printf "%2d examples compiled successfully in serial    " $nSerBld
if [[ $nSerErr > 0 ]]
then
	printf ", %2d failed" $nSerErr
fi
printf "\n%2d examples compiled successfully in parallel" $nParBld
if [[ $nParErr > 0 ]]
then
	printf ", %2d failed" $nParErr
fi
printf "\n%2d examples executed successfully            " $nParRun
if [[ $nRunErr > 0 ]]
then
	printf ", %2d failed" $nRunErr
fi
echo
cd ${root_dir}/test

AllERR=$(echo "$nSerErr+$nParErr+$nRunErr" | bc -l)
if [[ $AllERR > 0 ]]
then
	echo
	echo -e "${BRED}Build error(s) detected: ${AllERR} tests failed.${WHT}"
	exit 1
fi

echo
echo -e "${BGRN}Build tests successfully completed.${WHT}"
