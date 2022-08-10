compile:
	g++ spectralGap.cpp -o sG -I../eigen

clean:
	rm output/*
	rm output/Graphs/*

graph:
	cd output/ && gnuplot < ../spectralGap_commands.txt
	convert -delay 1 -loop 0 *.png Gap.gif
