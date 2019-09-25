echo Testing PH in epsilon range 0 - 3

for i in 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 
do
	echo "$i"
	../src/TDA_Cplusplus -i ../src/out.csv -e $i
done
