for i in {1..19}
do
	echo "Gen-$i"
	../GP -gen $i -gp best
	mv image.ps "Gen$i"_GPbest.ps
done
