
for i in $(find \( -name '*.cpp' -o -name '*.h' \) -a -not \( -path '*ThirdParty*' -o -path './.*' \) )
do
	echo $i
	# remove header
	sed '1,/*\// d' $i > tmp

	# concatenate 
	cat README.txt tmp > $i

	# remove temporary files
	rm tmp
done



