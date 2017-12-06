package util.landau;

import util.SpectraUtil;

public interface LandauLevelModel {

	public void fit(SpectraUtil.LandauLevelPeak[] peaks);
	public void fit(SpectraUtil.LandauLevelPeak[][] peaks);
	public double energy(double n, double B);
	public void setParam(double[] param);

}
