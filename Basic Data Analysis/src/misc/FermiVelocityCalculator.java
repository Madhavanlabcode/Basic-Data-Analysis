package misc;

public class FermiVelocityCalculator {

	static final double hbarJs = 1.054571726e-34;
	static final double hbarkgA2perS = hbarJs * 1e20;
	static final double eCoulomb = 1.60217657e-19;
	static final double hbarevs = hbarJs/eCoulomb;
	static final double e_ESU = 4.80320425e-10;
	static final double pi = Math.PI;
	/**
	 * This converts the momentum unit "sqrt(T)" to inverse Angstroms.
	 * One Tesla equals 1 kg s^-1 Coulomb^-1.
	 * 
	 * @param sqrtT
	 * @return
	 */
	public static double kInvA(double sqrtT)
	{
		return Math.sqrt(2*eCoulomb/hbarkgA2perS)*sqrtT;
	}
	public static double vectorPotential_teslaAngstroms(double kInvA){
		return (hbarkgA2perS/eCoulomb)*kInvA;
	}
	public static double fieldInInvA2(double bTesla){
		return bTesla*eCoulomb/hbarkgA2perS;
	}
	public static double kInvNM(double sqrtT)
	{
		return kInvA(sqrtT)*10;
	}
	/**
	 * This converts electronvolt Angstroms to meters per second.
	 * @param eva
	 * @return
	 */
	public static double mps(double eva)
	{
		double evm = eva/1e10;
		return evm/hbarevs;
	}
	
	/**
	 * evPerSqrtT is eV in the numerator and sqrtT in the denominator.
	 * @param evPerSqrtT
	 * @return
	 */
	public static double vf_mps(double evPerSqrtT)
	{
		double vf_eva = evPerSqrtT / kInvA(1);
		System.out.println(vf_eva);
		return mps(vf_eva);
	}
	
	
	public static void main (String[] args)
	{
		System.out.println(vf_mps(0.1/4));
//		System.out.println(mps(2.974));
//		System.out.println(kInvA(6));
//		int npts = 1024;
//		double emin = -500, emax = 100, ed0 = -120, ed = -200;
//		double[] v = ArrayOps.generateArrayInclBoth(emin, emax, npts);
//		double[] e = new double [npts];
//		for (int i = 0; i < npts; i++)
//			e[i] = v[i] + ed0*(1 - Math.sqrt(1 + (ed/ed0 - ed0/ed)*v[i]/ed0));
//		GraphDrawerCart.plotGraph(v, e);
		
	}
}
