
for i in 1 2 3 4 5 6 7 8 9 10; do
	for j in 1 2; do
		awk '{sub(/^[^[:space:]]+[[:space:]]*/, ""); print}' $1/profoundAutismBoth_rand_${i}_${j}.psam > $1/profoundAutismBoth_rand_${i}_${j}_noFamID.psam
	mv $1/profoundAutismBoth_rand_${i}_${j}.psam $1/profoundAutismBoth_rand_${i}_${j}_withFamilyID.psam
	mv $1/profoundAutismBoth_rand_${i}_${j}_noFamID.psam $1/profoundAutismBoth_rand_${i}_${j}.psam
	done
done
