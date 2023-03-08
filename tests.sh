./gaussfit -m -b -1 -s 5 -g A=10,sigma=2,mu=-3 -g ,A=10,sigma=2,mu=3 -p -f test_data/test3.dat
./gaussfit -s 90 -b -1 -m -g A=0.6,sigma=2,mu=328 -g A=0.7,sigma=7,mu=377 -g A=0.1,sigma=69,mu=383 -g A=0.04,sigma=15,mu=920 -g A=0.02,sigma=20,mu=1030 -p -f  test_data/test1.dat
./gaussfit -s 1 -b 0 -m -g A=250,sigma=20,mu=405 -g A=500,sigma=30,mu=440 -g A=550,sigma=69,mu=487 -g A=0.04,sigma=15,mu=920 -p -f  test_data/test2.dat
./gaussfit -s 1 -b 0 -m -l A=90,sigma=10,mu=405 -l A=500,sigma=30,mu=440 -l A=550,sigma=69,mu=467 -p -f  test_data/test2.dat
#Automatic mode, for now it only works in simple cases:
./gaussfit -a tyf=0 -b -1 -s 5 -p -f test_data/test3.dat
./gaussfit -a rth=0.18,tyf=1,nt=25,faf=1 -b -1 -s 5 -p -f test_data/test1.dat
./gaussfit -a rth=0.18,tyf=1,nt=25,faf=1 -g A=0.04,sigma=15,mu=920 -g A=0.02,sigma=20,mu=1030 -b -1 -s 5 -p -f test_data/test1.dat
./gaussfit -s 1 -b 0 -a null -p -f  test_data/test2.dat
