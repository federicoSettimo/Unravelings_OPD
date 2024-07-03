roqj.o: roqj.h roqj.cpp
	g++ roqj.cpp -c -o roqj.o -std=c++20 -O3 -ffast-math -fno-math-errno

roqj_pop.o: roqj_pop.cpp roqj_pop.h
	g++ roqj_pop.cpp -c -o roqj_pop.o -std=c++20 -O3 -ffast-math -fno-math-errno

JC_APO_continuum: roqj.o roqj_pop.o JC_continuum/JC_APO_continuum.cpp
	g++ JC_continuum/JC_APO_continuum.cpp roqj.o roqj_pop.o -o JC_continuum/JC_APO_continuum.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./JC_continuum/JC_APO_continuum.x
	#mv *.txt JC_continuum/.
	#python3 JC_continuum/plot_JC_APO.py

QSD_JC_APO_continuum: JC_continuum/QSD_JC_APO_continuum.cpp
	g++ JC_continuum/QSD_JC_APO_continuum.cpp -o JC_continuum/QSD_JC_APO_continuum.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./JC_continuum/QSD_JC_APO_continuum.x
	mv *.txt JC_continuum/.
	#python3 JC_continuum/plot_JC_APO.py