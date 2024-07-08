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

JC_APO_zero_discord: JC_single_mode/JC_APO_zero_discord.cpp
	g++ JC_single_mode/JC_APO_zero_discord.cpp -o JC_single_mode/JC_APO_zero_discord.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./JC_single_mode/JC_APO_zero_discord.x
	mv *.txt JC_single_mode/.
	#python3 JC_single_mode/plot_JC_APO.py

JC_APO_ent: JC_single_mode/JC_APO_ent.cpp
	g++ JC_single_mode/JC_APO_ent.cpp -o JC_single_mode/JC_APO_ent.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./JC_single_mode/JC_APO_ent.x
	mv *.txt JC_single_mode/.
	#python3 JC_single_mode/plot_JC_APO.py

sigma_xz: roqj.o roqj_pop.o sigma_xz/sigma_x_sigma_z.cpp
	g++ sigma_xz/sigma_x_sigma_z.cpp roqj.o roqj_pop.o -o sigma_xz/sigma_x_sigma_z.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./sigma_xz/sigma_x_sigma_z.x 0 p
	mv analytic.txt analytic_0_p.txt
	./sigma_xz/sigma_x_sigma_z.x 0 m
	mv analytic.txt analytic_0_m.txt
	./sigma_xz/sigma_x_sigma_z.x x p
	mv analytic.txt analytic_x_p.txt
	./sigma_xz/sigma_x_sigma_z.x x m
	mv analytic.txt analytic_x_m.txt
	./sigma_xz/sigma_x_sigma_z.x y p
	mv analytic.txt analytic_y_p.txt
	./sigma_xz/sigma_x_sigma_z.x y m
	mv analytic.txt analytic_y_m.txt
	mv *.txt sigma_xz/.

sigma_xz_Z: roqj.o roqj_pop.o sigma_xz/sigma_x_sigma_z_Z.cpp
	g++ sigma_xz/sigma_x_sigma_z_Z.cpp roqj.o roqj_pop.o -o sigma_xz/sigma_x_sigma_z_Z.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./sigma_xz/sigma_x_sigma_z_Z.x p
	mv analytic.txt analytic_z_p.txt
	./sigma_xz/sigma_x_sigma_z_Z.x m
	mv analytic.txt analytic_z_m.txt
	mv *.txt sigma_xz/.
	#python3 JC_single_mode/plot_JC_APO.py