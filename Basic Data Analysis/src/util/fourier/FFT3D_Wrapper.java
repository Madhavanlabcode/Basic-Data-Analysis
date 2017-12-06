package util.fourier;

import util.Complex;
import util.FieldOps;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_3D;

/**
 * This class provides an interface to the jtransforms methods for doing 3D transforms on topomaps.
 * @author madhavanlab2011
 *
 */
public class FFT3D_Wrapper {

	public double[][][] data;
	double[] a;
	public double[][][][] output; //This data is, of course, complex. [nlayers][nx][ny][2].
	DoubleFFT_3D jtFFT = null;
	
	int nlayers, nx, ny;
	int layerStride, xStride;
	
	public FFT3D_Wrapper(double[][][] data)
	{
		this.data = data;
		nlayers = data.length;
		nx = data[0].length;
		ny = data[0][0].length;
		
		layerStride = nx*ny*2;
		xStride = ny*2;
		
		jtFFT = new DoubleFFT_3D(nlayers, nx, ny);
	}
	
	public void doFFT()
	{
		a = new double [nlayers*nx*ny*2];
		for (int k = 0; k < nlayers; k++)
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
				{
					a[k*layerStride + i*xStride + 2*j] = data[k][i][j];
					a[k*layerStride + i*xStride + 2*j + 1] = 0;
				}
		jtFFT.complexForward(a);
		
		output = new double [nlayers][nx][ny][2];
		for (int k = 0; k < nlayers; k++)
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
				{
					output[k][i][j][0] = a[k*layerStride + i*xStride + 2*j];
					output[k][i][j][1] = a[k*layerStride + i*xStride + 2*j + 1];
				}
		a = null;
	}
	
	public double[][][] getFFTMag()
	{
		a = new double [nlayers*nx*ny*2];
		for (int k = 0; k < nlayers; k++)
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
				{
					a[k*layerStride + i*xStride + 2*j] = data[k][i][j];
					a[k*layerStride + i*xStride + 2*j + 1] = 0;
				}
		jtFFT.complexForward(a);
		
		double[][][] ans = new double [nlayers][nx][ny];
		for (int k = 0; k < nlayers; k++)
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
				{
					ans[k][i][j] = Complex.mag(a[k*layerStride + i*xStride + 2*j], a[k*layerStride + i*xStride + 2*j + 1]);
				}
		a = null;
		return ans;
	}
	
	//now, we put the COMPLEX result in output, and the real part into data
	public void doIFFT()
	{
		a = new double [nlayers*nx*ny*2];
		for (int k = 0; k < nlayers; k++)
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
				{
					a[k*layerStride + i*xStride + 2*j] = output[k][i][j][0];
					a[k*layerStride + i*xStride + 2*j + 1] = output[k][i][j][1];
				}
		jtFFT.complexInverse(a, true);
		
		if (output == null) output = new double [nlayers][nx][ny][2];
		for (int k = 0; k < nlayers; k++)
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
				{
					output[k][i][j][0] = a[k*layerStride + i*xStride + 2*j];
					output[k][i][j][1] = a[k*layerStride + i*xStride + 2*j + 1];
					data[k][i][j] = output[k][i][j][0];
				}
		a = null;
		
	}
	
	public static double[][][] getFFTMagCentXY(double[][][] data, boolean log)
	{
		FFT3D_Wrapper fft = new FFT3D_Wrapper(data);
		double[][][] ans = fft.getFFTMag();
		double[][] temp = new double [fft.nx][fft.ny];
		for (int i = 0; i < ans.length; i++)
		{
			FieldOps.shift(ans[i], temp);
			for (int j = 0; j < fft.nx; j++)
				for (int k = 0; k < fft.ny; k++)
					ans[i][j][k] = temp[j][k];
			
			if (log) FieldOps.log(ans[i]);
		}
		return ans;
	}
}
