compile:
	g++ -O3 --std=c++11 assignment-2018.c -o assignment
format:
	@echo "'1sf' timestamp total posX posY posZ velX velY velZ mass ... repeat from posX to add more body"
	@echo "ex1: 0.01  100.0   0 0 0 1.0   0   0 1.0"
	@echo "ex2: 0.01  100.0   0 0 0 1.0   0   0 1.0 0 1.0 0 1.0 0   0 1.0"
	@echo "ex3: 0.01  100.0 3.0 0 0   0 1.0   0 0.4 0   0 0   0 0   0 0.2 2.0 0 0 0 0 0 1.0"
run:
	./assignment 0.01  100.0 3.0 0 0   0 1.0   0 0.4 0   0 0   0 0   0 0.2 2.0 0 0 0 0 0 1.0
clean:
	rm *.vtp
	rm result.pvd
