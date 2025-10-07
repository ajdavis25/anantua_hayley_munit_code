for modelToDo in $(seq 0 1 9); do
	#10 grmhd models total.
	for kappaToDo in $(seq 0 1 2); do
		#4 edfs total.
		for sourceToDo in $(seq 0 1 1); do
			#2 sources total.
			echo ${modelToDo} ${kappaToDo} ${sourceToDo}
			sbatch fit_Munits_kappa_v3.sbatch ${modelToDo} ${kappaToDo} ${sourceToDo}
			sleep 0.5
		done
	done
done
