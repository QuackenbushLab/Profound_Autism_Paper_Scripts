
for i in 1 2 3 4 5 6 7 8 9 10; do
	for j in 1 2; do
		plink2 --bfile $3 --keep $1/profoundAutismBoth_rand_samps_${i}_${j}.csv --make-pgen --out $2/profoundAutismBoth_rand_${i}_${j}
	done
done
