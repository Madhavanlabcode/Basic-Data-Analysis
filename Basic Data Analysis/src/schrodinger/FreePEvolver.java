package schrodinger;

import util.FieldOps;
import util.fourier.FFT2D;

//This class determines the wave function at time t, assuming there is no potential.
//2D
public class FreePEvolver {

	//hbar and m are inputs. The unit of length is taken as the lattice spacing a
	//we define in natural units all quantities.
	double hbar, m, L;
	
	int N;
	
	double E0, p0, v0, t0, tau;

	//The particle will completely reproduce itself after a time equal to 
	//(4/pi)*tau
	double T;
	
	boolean stepConst = true;
	double dt;

	double[][][] psi;
	double[][][] psiHat0;
	double[][][] psiHat_t;
	double[][][] phases; //This applies if there is a constant time step dt.

	FFT2D fft; //this actually performs the necessary operations.
	
	//We assume it is square. dt is entered in natural units: that is, the user
	//is really entering dt/tau. It follows that for particles confined in boxes
	//of size ~L, the condition for short step will always be (input dt << 1).
	public FreePEvolver(double[][][] psi0, double dt, double hbar, double m)
	{
		N = psi0.length;
		
		//The physical constants. tau is given by mL^2/hbar;
		L = N; this.hbar = hbar; this.m = m;
		p0 = hbar/L; v0 = p0/m; E0 = hbar*hbar/(m*L*L);
		tau = hbar/E0;
		t0 = tau/(L*L); //t0 is the "actual" natural unit of time.
		T = 4*tau/Math.PI;
		this.dt = dt*tau;
		stepConst = true;
		
		//tau is given by mL^2/hbar;
		
		psi = new double[N][N][2];
		psiHat0 = new double [N][N][2];
		psiHat_t = new double [N][N][2];
		
		phases = new double [N][N][2];
		FieldOps.copy(psi0, psi);
		
		setPhases();
		
		fft = new FFT2D(psi);
		fft.doFFT();
		FieldOps.copy(fft.fHat, psiHat0);
		FieldOps.copy(psiHat0, psiHat_t);
	}
	
	void setPhases()
	{
		//This sets the phase factors which advance the wave function by dt.
		//The wavevector is pi/L (nx, ny)
		//the momentum is hbar*pi/L (nx, ny) = p0*pi (nx, ny)
		//the energy is p0^2 pi^2/2*m (nx, ny) = pi^2/2 E0 (nx^2 + ny^2)
		//The E_over_hbar is pi^2/2 * (1/tau) * (nx^2 + ny^2)
		//the largest E_over_hbar value will be ~ N^2 (1/tau)
		//dt << tau clearly insufficient if the higher modes are occupied.
		//in that case to fully capture the motion we should have dt << tau/N^2.
		//this however correspons to a small delta_x of order 1 instead of N.

		int nx, ny;
		double E_hbar;
		double pi2o2 = Math.PI*Math.PI/2;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				nx = i < N/2 ? i : i - N;
				ny = j < N/2 ? j : j - N;
				E_hbar = pi2o2 * (nx*nx + ny*ny)/tau;
				phases[i][j][0] = Math.cos(-E_hbar*dt);
				phases[i][j][1] = Math.sin(-E_hbar*dt);
			}
	}
	
	//advances psi by dt:
	void incrementPsiHat()
	{
		double[] temp = new double[2];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				product(psiHat_t[i][j], phases[i][j], psiHat_t[i][j], temp);
		fft.fHat = psiHat_t;
	}
	
	//This sets the value of psiHat at an arbitrary time. Note that this will overwrite
	//the incrementing and so should not be used simultaneuosly.
	void setPsiHat(double t)
	{
		double[] phase = new double[2];
		double[] temp = new double[2];
		int nx, ny;
		double E_hbar;
		double pi2o2 = Math.PI*Math.PI/2;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				nx = i < N/2 ? i : i - N;
				ny = j < N/2 ? j : j - N;
				E_hbar = pi2o2 * (nx*nx + ny*ny)/tau;
				phase[0] = Math.cos(-E_hbar*t);
				phase[1] = Math.sin(-E_hbar*t);
				product(psiHat0[i][j], phase, psiHat_t[i][j], temp);
			}
		
	}

	void obtainPsi()
	{
		fft.doIFFT();
	}
	
	//This returns psiHat centered at N/2, N/2
	public void getPsiHatCentered(double[][][] target)
	{
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				target[i][j][0] = psiHat_t[(i+N/2)%N][(j+N/2)%N][0];
				target[i][j][1] = psiHat_t[(i+N/2)%N][(j+N/2)%N][1];
			}
	}
	//This returns psi centered at N/2, N/2
	public void getPsiCentered(double[][][] target)
	{
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				target[i][j][0] = psi[(i+N/2)%N][(j+N/2)%N][0];
				target[i][j][1] = psi[(i+N/2)%N][(j+N/2)%N][1];
			}
	}
	
	//in case the target is either a or b.
	static void product(double[] a, double[] b, double[] target, double[] temp)
	{
		temp[0] = a[0]*b[0] - a[1]*b[1];
		temp[1] = a[0]*b[1] + a[1]*b[0];
		target[0] = temp[0];
		target[1] = temp[1];
	}
}
