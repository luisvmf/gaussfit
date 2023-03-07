# gaussfit

Decompose spectra in gaussians. See tests.sh for examples.

		Usage:
		 gaussfit [-p] [-m] [-a] [-g A=XX,sigma=XX,mu=XX] [-g A=XX,sigma=XX,mu=XX] [-f file]

		   -w X  File save mode.
				     X=0-> Don't save results (default).
				     X=1-> Saves final fit parameters
				           (filename.gaussfit).
				     X=2-> Saves final fit parameters
				           and spectrum with boxcar and
				           baselined removed (filename.gaussfit
				           and filename.proc).
		   -b X  Baseline polynomial order (X is an integer).
		   -s X  Boxcar smoothing (X is an integer).
		   -p    Plot graph and gaussians with gnuplot.
				 Gnuplot must be installed and on path
				 for this option to be used.
		   -m    Manualy specify gaussians/Lorentzians initial
				 guess with option -g or -l. Must not be used
				 together with option -a.
		   -a    Automaticaly set gaussians initial guess.
		   -g A=XX,sigma=XX,mu=XX
				 Set gaussians initial guess, where XX are
				 float point values. Options -g or -l must
				 be used if option -m is used. 
		   -l A=XX,sigma=XX,mu=XX
				 Set Lorentzian initial guess, where XX are
				 float point values. Options -g or -l must
				 be used if option -m is used. 
		   -f file
				 Set file containing data to
				 be fitted. Must be used.

