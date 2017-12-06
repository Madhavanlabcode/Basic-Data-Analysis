package util.calc;

public interface Calculation2D_2DP {

	public double[][] realField();
	public double[][][] complexField();
	public double[][] outputField(int x, int y);
	
	public String getInfo(int x, int y);
	public String getInfo(double x, double y);
	
	public double[] getMinMax();
	public double[] getOutputMinMax();
	public void redoCalc();
	
	public void resize(double xmin, double xmax, double ymin, double ymax);
	public double[] getFieldBounds();
	
	public void switchDisplayField();
	public void save(String dir);

}
