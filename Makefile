compile:
	g++ spectralGap.cpp -o sG -I../eigen

clean:
	rm output/*
	rm output/Graphs/*

graph:
	cd output/ && gnuplot < ../spectralGap_commands.txt && rm Graphs/Graphs.png

reverse:
	cp output/Graphs/Gap.gif ./Gap.gif && gifsicle --colors=255 Gap.gif -o Gap-simplified.gif && gifsicle --unoptimize Gap-simplified.gif Gap-simplified.gif "#-1-0" > reversed.gif
	
video:
	cd output && ffmpeg -y -i ./Graphs/*.png -c:v libx264 -framerate 1/5 -r 30 -pix_fmt yuv420p out.mp4 
