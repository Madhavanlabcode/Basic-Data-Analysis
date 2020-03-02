package util;

import java.util.ArrayList;

import org.apache.commons.math3.analysis.interpolation.BicubicSplineInterpolatingFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import drawing.LayerViewer;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;
import flanagan.analysis.Regression;
import util.fileops.Layer;
import util.fourier.FFTOps;
import util.geom.AtomicCoordinatesSet;
import util.geom.Distance;
import util.matrix.Matrix;
import main.DataManip;

public class FieldOps {

//	public static double[][][] gradient(double[][] field)
//	{
//		int nx = field.length, ny = field[0].length;
//		double[][][] grad = new double [2][nx][ny]; //0 x, 1 y.
//		
//		int i = 0, j = 0;
//		//edges
//		for (i = 0; i < ny; i++)
//		{
//			grad[0][0][i] = field[1][i] - field[0][i];
//			grad[0][nx - 1][i] = field[nx - 1][i] - field[nx - 2][i];
//			if (i > 0 && i < ny-1)
//			{
//				grad[1][0][i] = field[0][i+1] - field[0][i-1];
//				grad[1][nx-1][i] = field[nx-1][i+1] - field[nx-1][i-1];
//			}
//		}
//		for (i = 0; i < nx; i++)
//		{
//			grad[1][i][0] = field[i][1] - field[i][0];
//			grad[1][i][ny - 1] = field[i][ny - 1] - field[i][ny - 2];
//			if (i > 0 && i < ny-1)
//			{
//				grad[0][i][0] = field[i+1][0] - field[i-1][0];
//				grad[0][ny-1][i] = field[i+1][ny-1] - field[i-1][ny-1];
//			}
//		}
//		
//		for (i = 1; i < nx - 1; i++)
//			for (j = 1; j < ny - 1; j++)
//			{
//				grad[0][i][j] = (field[i+1][j] - field[i-1][j])/2;
//				grad[1][i][j] = (field[i][j+1] - field[i][j-1])/2;
//			}
//		
//		return grad;
//	}
	
	public static double[][][] gradient(double[][] field)
	{
		int nx = field.length, ny = field[0].length;
		double[][][] grad = new double [2][nx][ny]; //0 x, 1 y.
		
		//First do the x derivative. The x derivative is ok except on the left and right edges.
		for (int i = 1; i < nx-1; i++)
			for (int j = 0; j < ny; j++)
				grad[0][i][j] = (field[i+1][j] - field[i-1][j])/2;
		
		for (int j = 0; j < ny; j++)
		{
			grad[0][0][j] = (field[1][j] - field[0][j]);
			grad[0][nx-1][j] = (field[nx-1][j] - field[nx-2][j]);
		}
		
		//Now the y-derivative, OK except at top and bottom
		for (int i = 0; i < nx; i++)
			for (int j = 1; j < ny-1; j++)
				grad[1][i][j] = (field[i][j+1]-field[i][j-1])/2;
		
		for (int i = 0; i < nx; i++)
		{
			grad[1][i][0] = (field[i][1]-field[i][0]);
			grad[1][i][ny-1] = (field[i][ny-1] - field[i][ny-2]);
		}
		
		return grad;
	}
	public static double[][] directionalDerivative2D(double[][] field, double[] direction)
	{
		int nx = field.length, ny = field[0].length;
		double[][][] grad = gradient(field);
		
		double[][] ans = new double [nx][ny];
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				ans[i][j] = direction[0]*grad[0][i][j] + direction[1]*grad[1][i][j];
			}
			
		return ans;
	}
	public static double[][] directionalDerivative2D(double[] direction, double[][][] grad)
	{
		int nx = grad[0].length, ny = grad[0][0].length;
//		double[][][] grad = gradient(field);
		
		double[][] ans = new double [nx][ny];
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				ans[i][j] = direction[0]*grad[0][i][j] + direction[1]*grad[1][i][j];
			}
			
		return ans;
	}
	public static double[][][] getDirectional2ndDerivativeStuff(double[][] field, double[] ahat, double[] bhat)
	{
		int nx = field.length, ny = field[0].length;
		double[][][] grad = gradient(field);
		double[][][] dirgrad = new double[][][] {directionalDerivative2D(ahat, grad), directionalDerivative2D(bhat, grad)};
		double[][][] temp = gradient(dirgrad[0]);
		double[][] faa = directionalDerivative2D(ahat, temp);
		double[][] fab = directionalDerivative2D(bhat, temp);
		temp = gradient(dirgrad[1]);
		double[][] fba = directionalDerivative2D(ahat, temp);
		double[][] fbb = directionalDerivative2D(bhat, temp);
		
		
		double[][] e12 = new double [nx][ny];
		double[][] w12 = new double [nx][ny];
		double[][] tr_e = new double [nx][ny];
		double[][] anis = new double [nx][ny];
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				e12[i][j] = (fab[i][j] + fba[i][j])/2;
				w12[i][j] = (fab[i][j] - fba[i][j])/2;
				tr_e[i][j] = (faa[i][j] + fbb[i][j])/2;
				anis[i][j] = (faa[i][j] - fbb[i][j])/2;
			}
		
		return new double[][][] {faa, fab, fba, fbb, e12, w12, tr_e, anis};
	}
	
	public static double[][] getSecondDirectionalDerivative2D(double[][] field, double[] direction)
	{
		return directionalDerivative2D(directionalDerivative2D(field, direction), direction);
	}
	
	/**
	 * This calculates the curvature thing described in  Rev. Sci. Instrum. 82, 043712 (2011); http://dx.doi.org/10.1063/1.3585113 
	 * where the parameter is given in units of the maximum value of the squared magnitude of the gradient. 
	 * Large parameter is like the laplacian, small parameter is like one 
	 * The parameter 
	 * @param data
	 * @param parameter
	 * @return
	 */
	public static double[][] getCurvatureThing(double[][] f, double parameter)
	{
		double[][][] grad = gradient(f);
		double[][] dfdx = grad[0];
		double[][] dfdy = grad[1];
		double[][][] grad_x = gradient(dfdx);
		
		
		
		double[][] d2fdx2 = grad_x[0];
		double[][] d2fdxdy = grad_x[1];
		
		double[][] d2fdy2 = gradient(dfdy)[1];
		
		double[][] ans = new double [f.length][f[0].length];
		
		double gradMagMax = FieldOps.magMax(dfdx, dfdy);
		
		double c0 = parameter*gradMagMax*gradMagMax;
//		LayerViewer.show(Layer.getFreeLayer(dfdx), 1024);
//		LayerViewer.show(Layer.getFreeLayer(dfdy), 1024);
//		LayerViewer.show(Layer.getFreeLayer(d2fdx2), 1024);
//		LayerViewer.show(Layer.getFreeLayer(d2fdxdy), 1024);
//		LayerViewer.show(Layer.getFreeLayer(d2fdy2), 1024);
//		LayerViewer.show(Layer.getFreeLayer(FieldOps.magnitude(fftz)), 1024);
//		LayerViewer.show(Layer.getFreeLayer(FieldOps.magnitude(fftz)), 1024);
//		LayerViewer.show(Layer.getFreeLayer(FieldOps.magnitude(fftz)), 1024);
//		LayerViewer.show(Layer.getFreeLayer(FieldOps.magnitude(fftz)), 1024);

		for (int i = 0; i < f.length; i++)
			for (int j = 0; j < f[0].length; j++)
			{
//				ans[i][j] = ( (c0 + Math.pow(dfdx[i][j], 2))*d2fdy2[i][j] - (2*dfdx[i][j]*dfdy[i][j]*d2fdxdy[i][j]) + (c0+Math.pow(dfdy[i][j],2)*d2fdx2[i][j]))
				ans[i][j] = (c0*(d2fdy2[i][j] + d2fdx2[i][j]) + d2fdy2[i][j]*Math.pow(dfdx[i][j], 2) + d2fdx2[i][j]*Math.pow(dfdy[i][j], 2) - (2*dfdx[i][j]*dfdy[i][j]*d2fdxdy[i][j]))
				/
						Math.pow(c0 + Math.pow(dfdx[i][j], 2) + Math.pow(dfdy[i][j], 2), 1.5);
			}
		return ans;
	}
	//takes the laplacian in the form laplacian = inverse fourier transorm of -k^2(fHat(k))
	public static double[][] laplaceFourier(double[][] field)
	{
		
		int nx = field.length, ny = field[0].length;
		double[][][] fftz = new double [nx][ny][2];
		FFTOps.putFFT(field, fftz, true);
//		LayerViewer.show(Layer.getFreeLayer(FieldOps.magnitude(fftz)), 1024);
		double kx, ky;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				kx = i-nx/2;
				ky = j-ny/2;
				if (kx != 0 || ky != 0){
					fftz[i][j][0] *= -(kx*kx + ky*ky);
					fftz[i][j][1] *= -(kx*kx + ky*ky);
				}
			}
//		LayerViewer.show(Layer.getFreeLayer(FieldOps.magnitude(fftz)), 1024);
		double[][][] ifft = new double [nx][ny][2];
		FFTOps.putIFFT(fftz, ifft, true);
		double[][] ans = new double [nx][ny];
		FieldOps.getIndex(ifft, ans, 0);
		return ans; 
	}
	public static double[][] inverseLaplaceFourier(double[][] field)
	{
		
		int nx = field.length, ny = field[0].length;
		double[][][] fftz = new double [nx][ny][2];
		FFTOps.putFFT(field, fftz, true);
		double kx, ky;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				kx = i-nx/2;
				ky = j-ny/2;
				if (kx != 0 || ky != 0){
					fftz[i][j][0] /= -(kx*kx + ky*ky);
					fftz[i][j][1] /= -(kx*kx + ky*ky);
				}
			}
		double[][][] ifft = new double [nx][ny][2];
		FFTOps.putIFFT(fftz, ifft, true);
		double[][] ans = new double [nx][ny];
		FieldOps.getIndex(ifft, ans, 0);
		return ans; 
	}
	public static double[][] inverseLaplaceReal(double[][] field)
	{
		
		int nx = field.length, ny = field[0].length;
		double[][] ans = new double [nx][ny];
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				for (int p = 0; p < nx; p++)
					for (int q = 0; q < ny; q++)
					{
						double r = Distance.distance(i-p, j-q);
						if (r != 0)
							ans[i][j] += field[p][q]*Math.log(r);
					}
		return ans; 
	}
	public static double[][] laplace(double[][] field)
	{
		
		int nx = field.length, ny = field[0].length;
		double[][] ans = new double [nx][ny];
		for (int i = 1; i < nx-1; i++)
			for (int j = 1; j < ny-1; j++)
			{
				ans[i][j] = field[i+1][j]+field[i-1][j] + field[i][j+1] + field[i][j-1] - 4*field[i][j];
			}
		for (int j = 1; j < ny-1; j++)
		{
			ans[0][j] = field[0][j+1] + field[0][j-1] - 2*field[0][j];
			ans[nx-1][j] = field[nx-1][j+1] + field[nx-1][j-1] - 2*field[nx-1][j];
		}
		for (int i = 1; i < nx-1; i++)
		{
			ans[i][0] = field[i+1][0]+field[i-1][0] - 2*field[i][0];
			ans[i][ny-1] = field[i+1][ny-1] + field[i-1][ny-1] - 2*field[i][ny-1];
		}
		//ignore the edges for now
		return ans; 
	}
	public static double[][][] laplace(double[][][] field)
	{
		
		int nx = field.length, ny = field[0].length, nk = field[0][0].length;
		double[][][] ans = new double [nx][ny][nk];
		for (int i = 1; i < nx-1; i++)
			for (int j = 1; j < ny-1; j++)
				for (int k = 1; k < nk-1; k++)
				{
					ans[i][j][k] = field[i+1][j][k]+field[i-1][j][k] + field[i][j+1][k] + field[i][j-1][k] + field[i][j][k-1] + field[i][j][k+1] - 6*field[i][j][k];
				}
		/**
		 * This must now handle the planes i = 0 and i = nx-1
		 */
		for (int j = 1; j < ny-1; j++)
			for (int k = 1; k < nk-1; k++)
			{
				ans[0][j][k] = field[0][j+1][k] + field[0][j-1][k] + field[0][j][k+1] + field[0][j][k-1] - 4*field[0][j][k];
				ans[0][j][k] += field[2][j][k] + field[0][j][k] - 2*field[1][j][k];
				ans[nx-1][j][k] = field[nx-1][j+1][k] + field[nx-1][j-1][k] + field[nx-1][j][k+1] + field[nx-1][j][k-1] - 4*field[nx-1][j][k];
				ans[nx-1][j][k] += field[nx-1][j][k] + field[nx-3][j][k] - 2*field[nx-2][j][k];
			}
		/**
		 * The planes j=0 and j=ny-1
		 */
		for (int i = 1; i < nx-1; i++)
			for (int k = 1; k < nk-1; k++)
			{
				ans[i][0][k] = field[i+1][0][k] + field[i-1][0][k] + field[i][0][k+1] + field[i][0][k-1] - 4*field[i][0][k];
				ans[i][0][k] += field[i][2][k] + field[i][0][k] - 2*field[i][1][k];
				ans[i][ny-1][k] = field[i+1][ny-1][k] + field[i-1][ny-1][k] + field[i][ny-1][k+1] + field[i][ny-1][k-1] - 4*field[i][ny-1][k];
				ans[i][ny-1][k] += field[i][ny-3][k] + field[i][ny-1][k] - 2*field[i][ny-2][k];
			}
		/**
		 * The planes k=0 and k = nk-1
		 *
		 */
		for (int i = 1; i < nx-1; i++)
			for (int j = 1; j < ny-1; j++)
			{
				ans[i][j][0] = field[i+1][j][0] + field[i-1][j][0] + field[i][j-1][0] + field[i][j+1][0] - 4*field[i][j][0];
				ans[i][j][0] += field[i][j][2] + field[i][j][0] - 2*field[i][j][1];
				ans[i][j][nk-1] = field[i+1][j][nk-1] + field[i-1][j][nk-1] + field[i][j-1][nk-1] + field[i][j+1][nk-1] - 4*field[i][j][nk-1];
				ans[i][j][nk-1] += field[i][j][nk-3] + field[i][j][nk-1] - 2*field[i][j][nk-2];
			}
		

		return ans; 
	}
	public static void putLaplaceFourier(double[][] field, double[][] ans)
	{
		
		int nx = field.length, ny = field[0].length;
		double[][][] fftz = new double [nx][ny][2];
		FFTOps.putFFT(field, fftz, true);
		double kx, ky;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				kx = i-nx/2;
				ky = j-ny/2;
				fftz[i][j][0] *= -(kx*kx + ky*ky);
				fftz[i][j][1] *= -(kx*kx + ky*ky);
			}
		double[][][] ifft = new double [nx][ny][2];
		FFTOps.putIFFT(fftz, ifft, true);
		FieldOps.getIndex(ifft, ans, 0);
	}
	public static double[][][] gradientNM2(double[][] field)
	{
		int nx = field.length, ny = field[0].length;
		double[][][] grad = new double [nx][ny][2]; //0 x, 1 y.
		
		int i = 0, j = 0;
		//edges
		for (i = 0; i < ny; i++)
		{
			grad[0][i][0] = field[1][i] - field[0][i];
			grad[nx - 1][i][0] = field[nx - 1][i] - field[nx - 2][i];
		}
		for (i = 0; i < nx; i++)
		{
			grad[i][0][1] = field[i][1] - field[i][0];
			grad[i][ny - 1][1] = field[i][ny - 1] - field[i][ny - 2];
		}
		
		for (i = 1; i < nx - 1; i++)
			for (j = 1; j < ny - 1; j++)
			{
				grad[i][j][0] = (field[i+1][j] - field[i-1][j])/2;
				grad[i][j][1] = (field[i][j+1] - field[i][j-1])/2;
			}
		
		return grad;
	}
	public static double[][] divergence(double[][][] vectorNM2)
	{
		int nx = vectorNM2.length, ny = vectorNM2[0].length;
		double[][] xComponent = FieldOps.getIndex(vectorNM2, 0);
		double[][] yComponent = FieldOps.getIndex(vectorNM2, 1);
		
		double[][][] gradX = FieldOps.gradient(xComponent); //0 x, 1 y.
		double[][][] gradY = FieldOps.gradient(yComponent);
		double[][] div = FieldOps.add(gradX[0], gradY[1]);
		return div;
	}
	public static double[][] divergence_2NM(double[][][] array)
	{
//		int nx = vectorNM2.length, ny = vectorNM2[0].length;
		double[][] xComponent = array[0];
		double[][] yComponent = array[1];
		
		double[][][] gradX = FieldOps.gradient(xComponent); //0 x, 1 y.
		double[][][] gradY = FieldOps.gradient(yComponent);
		double[][] div = FieldOps.add(gradX[0], gradY[1]);
		return div;
	}
	
	//This returns a scalar field, being the z-component of the curl of the two-compoenent 2-d vector field f.
	//we have f[i][j][0 or 1];
	public static double[][] curl(double[][][] f)
	{
		int nx = f.length, ny = f[0].length;
		double[][] curl = new double [nx][ny]; //0 x, 1 y.
		
		int i = 0, j = 0;
		//edges
		double dfydx, dfxdy;
		for (i = 1; i < ny-1; i++)
		{
			dfydx = f[1][i][1] - f[0][i][1];
			dfxdy = (f[0][i+1][0] - f[0][i-1][0])/2;
			curl[0][i] = dfydx - dfxdy;

			dfydx = f[nx-1][i][1] - f[nx-2][i][1];
			dfxdy = (f[nx-1][i+1][0] - f[nx-1][i-1][0])/2;
			curl[nx-1][i] = dfydx - dfxdy;
		}
		for (i = 1; i < nx-1; i++)
		{
			dfxdy = f[i][1][0] - f[i][0][0];
			dfydx = (f[i+1][0][1] - f[i-1][0][1])/2;
			curl[i][0] = dfydx - dfxdy;
			
			dfxdy = f[i][ny-1][0] - f[i][ny-2][0];
			dfydx = (f[i+1][ny-1][1] - f[i-1][ny-1][1])/2;
			curl[i][ny-1] = dfydx - dfxdy;
		}
		
		curl[0][0] = (f[1][0][1] - f[0][0][1]) - (f[0][1][0] - f[0][0][0]);
		curl[nx-1][0] = (f[nx-1][0][1] - f[nx-2][0][1]) - (f[nx-1][1][0] - f[nx-1][0][0]);
		curl[0][ny-1] = (f[1][ny-1][1] - f[0][ny-1][1]) - (f[0][ny-1][0] - f[0][ny-2][0]);
		curl[nx-1][ny-1] = (f[nx-1][ny-1][1] - f[nx-2][ny-1][1]) - (f[nx-1][ny-1][0] - f[nx-1][ny-2][0]);
		
		for (i = 1; i < nx - 1; i++)
			for (j = 1; j < ny - 1; j++)
			{
				dfxdy = (f[i][j+1][0] - f[i][j-1][0])/2;
				dfydx = (f[i+1][j][1] - f[i-1][j][1])/2;
				curl[i][j] = dfydx - dfxdy;
			}
		
		return curl;
	}
	public static double[][] curlSimple(double[][][] f)
	{
		int nx = f.length, ny = f[0].length;
		double[][] curl = new double [nx][ny]; //0 x, 1 y.
		double[][] dfydx = gradient(FieldOps.getIndex(f, 1))[0];
		double[][] dfxdy = gradient(FieldOps.getIndex(f, 0))[1];
		
		double[][] ans = FieldOps.minus(dfydx, dfxdy);
		return ans;
	}
	
	public static double[][][] getPseudofieldStuff(double[][] uaa, double[][] ubb, double[][] e12, double a1, double a2, double a3){
		int nx = uaa.length, ny = uaa[0].length;
		double[][][] A1 = new double [nx][ny][2];
		double[][][] A2 = new double [nx][ny][2];
		double[][] B1 = new double [nx][ny];
		double[][] B2 = new double [nx][ny];
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++){
				A1[i][j][0] = a1*uaa[i][j] + a2*ubb[i][j];
				A1[i][j][1] = a3*e12[i][j];
				A2[i][j][0] = a3*e12[i][j];
				A2[i][j][1] = a1*ubb[i][j] + a2*uaa[i][j];
			}
		
		B1 = curl(A1);
		B2 = curl(A2);
		return new double[][][] {getIndex(A1, 0), getIndex(A1, 1), getIndex(A2, 0), getIndex(A2, 1), B1, B2};
	}
	public static double[][][] getPseudofieldStuff_alternate(double[][] comp, double[][] anis, double[][] e12, double aC, double aU, double a3){
		int nx = comp.length, ny = comp[0].length;
		double[][][] A1 = new double [nx][ny][2];
		double[][][] A2 = new double [nx][ny][2];
		double[][] B1 = new double [nx][ny];
		double[][] B2 = new double [nx][ny];
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++){
				A1[i][j][0] = aC*comp[i][j] + aU*anis[i][j];
				A1[i][j][1] = a3*e12[i][j];
				A2[i][j][0] = a3*e12[i][j];
				A2[i][j][1] = aC*comp[i][j] - aU*anis[i][j];
			}
		
		B1 = curlSimple(A1);
		B2 = curlSimple(A2);
		return new double[][][] {getIndex(A1, 0), getIndex(A1, 1), getIndex(A2, 0), getIndex(A2, 1), B1, B2};
	}
	public static void abs(double[][] f)
	{
		for (int i = 0; i < f.length; i++)
			for (int j = 0; j < f.length; j++)
				f[i][j] = Math.abs(f[i][j]);
	}
	public static double[][] getAbs(double[][] f)
	{
		double[][] g = new double [f.length][f[0].length];
		for (int i = 0; i < f.length; i++)
			for (int j = 0; j < f.length; j++)
				g[i][j] = Math.abs(f[i][j]);
		return g;
	}
	public static void star(double[][][] source, double[][][] target)
	{
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
			{
				target[i][j][0] = source[i][j][0];
				target[i][j][1] = -source[i][j][1];
			}
	}
	public static void copy(double[][][] source, double[][][] target)
	{
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[i].length; j++)
				for (int k = 0; k < source[i][j].length; k++)
				{
					target[i][j][k] = source[i][j][k];
				}
	}
	public static double[][][] copy(double[][][] source)
	{
		double[][][] target = new double [source.length][source[0].length][source[0][0].length];
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[j].length; j++)
				for (int k = 0; k < source[i][j].length; k++)
				{
					target[i][j][k] = source[i][j][k];
				}
		return target;
	}
	public static double[][][][] copy(double[][][][] source)
	{
		double[][][][] target = new double [source.length][source[0].length][source[0][0].length][source[0][0][0].length];
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[i].length; j++)
				for (int k = 0; k < source[i][j].length; k++)
					for (int l = 0; l < source[i][j][k].length; l++)
				{
					target[i][j][k][l] = source[i][j][k][l];
				}
		return target;
	}
	public static double[] copy(double[] source)
	{
		double[] target = new double [source.length];
		for (int i = 0; i < source.length; i++)
				target[i]=source[i];
		return target;
	}
	public static void copy(double[][] source, double[][] target)
	{
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
			{
				target[i][j] = source[i][j];
				target[i][j] = source[i][j];
			}
	}
	public static void copy(double[][] source, double[][][] target, int target1stindex)
	{
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
			{
				target[target1stindex][i][j] = source[i][j];
				target[target1stindex][i][j] = source[i][j];
			}
	}
	public static double[][] copy(double[][][] source, int thirdIndex)
	{
		double[][] target = new double [source.length][source[0].length];
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
			{
				target[i][j] = source[i][j][thirdIndex];
			}
		return target;
	}
	public static boolean[][] copy(boolean[][] source)
	{
		boolean[][] target = new boolean [source.length][source[0].length];
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
			{
				target[i][j] = source[i][j];
			}
		return target;
	}
	public static void copy(int[][] source, double[][] target)
	{
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
			{
				target[i][j] = source[i][j];
				target[i][j] = source[i][j];
			}
	}
	public static double[][] truncate(double[][] source, int imin, int imax, int jmin, int jmax)
	{
		double[][] out = new double [imax-imin][jmax-jmin];
		for (int i = imin; i < imax; i++)
			for (int j = jmin; j < jmax; j++)
			{
				out[i-imin][j-jmin] = source[i][j];
			}
		return out;
	}
	public static double[][] augment(double[][] source, double[] ilabels, double[] jlabels)
	{
		double[][] out = new double [source.length+1][source[0].length+1];
		out[0][0] = Double.NaN;
		for (int i = 0; i < source.length; i++)
			out[i+1][0] = ilabels[i];
		for (int j = 0; j < source[0].length; j++)
			out[0][j+1] = jlabels[j];
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
			{
				out[i+1][j+1] = source[i][j];
				out[i+1][j+1] = source[i][j];
			}
		return out;
	}
	
	public static double[][] normalizeEachRow(double[][] in, boolean firstIndexRow)
	{
		double[][] out = new double [in.length][in[0].length];
		double max;
		if (firstIndexRow)
			for (int i = 0; i < in[0].length; i++)
			{
				max = ArrayOps.max(in, i, firstIndexRow);
				for (int j = 0; j < in.length; j++)
					out[j][i] = in[j][i]/max;
			}
		else
			for (int i = 0; i < in.length; i++)
			{
				max = ArrayOps.max(in, i, firstIndexRow);
				for (int j = 0; j < in[0].length; j++)
					out[i][j] = in[i][j]/max;
			}
		return out;
	}
	//normalizes by dividing by the sum of squared values, as in quantum mechanics.
	public static double[][] normalizeWaveFunction(double[][] in)
	{
		double[][] out = new double [in.length][in[0].length];
		double integral = 0;
			for (int i = 0; i < in[0].length; i++)
				for (int j = 0; j < in.length; j++)
					integral += in[i][j]*in[i][j];
			
			for (int i = 0; i < in.length; i++)
				for (int j = 0; j < in[0].length; j++)
					out[i][j] = in[i][j]/integral;
				
			return out;
	}
	public static void mask(double[][] source, double[][] target, int[][] region, int key)
	{
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
				if (region[i][j] == key) target[i][j] = source[i][j];
				else target[i][j] = 0;
	}
	public static void mask(double[][] source, double[][] target, int[][] region, int key, double val)
	{
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
				if (region[i][j] == key) target[i][j] = source[i][j];
				else target[i][j] = val;
	}
	public static double getMeanRegion(double[][] source, int[][] region, int key)
	{
		double m = 0; int n = 0;
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
				if (region[i][j] == key) {m += source[i][j]; n++;}
		return m/n;
	}
	public static double[][] copy(double[][] source)
	{
		double[][] target = new double [source.length][source[0].length];
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
			{
				target[i][j] = source[i][j];
				target[i][j] = source[i][j];
			}
		return target;
	}
	
	public static double[][] magnitude(double[][] x, double[][] y)
	{
		int N = x.length, M = x[0].length;
		double[][] mag = new double[N][M];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				mag[i][j] = Math.sqrt(x[i][j]*x[i][j] + y[i][j]*y[i][j]);
		return mag;
	}
	public static double[][] phase(double[][] x, double[][] y)
	{
		int N = x.length, M = x[0].length;
		double[][] phase = new double[N][M];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				phase[i][j] = atan(x[i][j], y[i][j]);
		return phase;
	}
	public static double[][] magnitude(double[][][] f)
	{
		int N = f.length, M = f[0].length;
		double[][] mag = new double[N][M];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				mag[i][j] = Math.sqrt(f[i][j][0]*f[i][j][0] + f[i][j][1]*f[i][j][1]);
		return mag;
	}
	
	public static void magnitude(double[][][] f, double[][] target)
	{
		int N = f.length, M = f[0].length;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				target[i][j] = Math.sqrt(f[i][j][0]*f[i][j][0] + f[i][j][1]*f[i][j][1]);
	}
	public static void getIndex(double[][][] f, double[][] target, int thirdindex)
	{
		int N = f.length, M = f[0].length;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				target[i][j] = f[i][j][thirdindex];
	}
	public static double[][] getIndex(double[][][] f, int thirdindex)
	{
		int N = f.length, M = f[0].length;
		double[][] target = new double [N][M];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				target[i][j] = f[i][j][thirdindex];
		return target;
	}
	public static void putIntoArray(double[][] ux, double[][] uy, double[][][] f)
	{
		int N = f.length, M = f[0].length;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++){
				f[i][j][0] = ux[i][j];
				f[i][j][1] = uy[i][j];
			}
	}
	public static double[][] phase(double[][][] f)
	{
		int N = f.length, M = f[0].length;
		double[][] phase = new double[N][M];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				phase[i][j] = atan(f[i][j][0], f[i][j][1]);
		return phase;
	}
	public static double[][] phaseCenteredZero(double[][][] f)
	{
		int N = f.length, M = f[0].length;
		double[][] phase = new double[N][M];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++){
				phase[i][j] = atan(f[i][j][0], f[i][j][1]);
				if (phase[i][j] > Math.PI) phase[i][j] -= 2*Math.PI;
			}
		return phase;
	}
	public static void phase(double[][][] f, double[][] target)
	{
		int N = f.length, M = f[0].length;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				target[i][j] = atan(f[i][j][0], f[i][j][1]);
	}
	
	public static void plusEquals(double[][] f, double delta)
	{
		int N = f.length, M = f[0].length;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				f[i][j] += delta;
	}
	public static void plusEquals(int[][] f, int delta)
	{
		int N = f.length, M = f[0].length;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				f[i][j] += delta;
	}
	public static void overEquals(double[][] f, double[][] x)
	{
		int N = f.length, M = f[0].length;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				f[i][j] /= x[i][j];
	}
	public static void minusEquals(double[][] f, double[][] x)
	{
		int N = f.length, M = f[0].length;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				f[i][j] -= x[i][j];
	}
	public static void minusEquals(double[][][] f, double[][][] x)
	{
		int N = f[0].length, M = f[0][0].length;
		for (int k = 0; k < f.length; k++)
			for (int i = 0; i < N; i++)
				for (int j = 0; j < M; j++)
					f[k][i][j] -= x[k][i][j];
	}
	public static void plusEquals(double[][][] f, double[][][] x)
	{
		int N = f[0].length, M = f[0][0].length;
		for (int k = 0; k < f.length; k++)
			for (int i = 0; i < N; i++)
				for (int j = 0; j < M; j++)
					f[k][i][j] += x[k][i][j];
	}
	public static void timesEquals(double[][] f, double x)
	{
		int N = f.length, M = f[0].length;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				f[i][j] *= x;
	}
	public static void minusEquals(double[][] f, double x)
	{
		int N = f.length, M = f[0].length;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				f[i][j] -= x;
	}
	public static void minusEquals(int[][] f, int x)
	{
		int N = f.length, M = f[0].length;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				f[i][j] -= x;
	}
	public static double[][] minus(double[][] f, double[][] x)
	{
		int N = f.length, M = f[0].length;
		double[][] a = new double [N][M];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				a[i][j] = f[i][j] - x[i][j];
		return a;
	}
	public static double[][] minus(double[][] f, double x)
	{
		int N = f.length, M = f[0].length;
		double[][] a = new double [N][M];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				a[i][j] = f[i][j] - x;
		return a;
	}
	public static void minus(double[][] f, double[][] x, double[][] target)
	{
		int N = f.length, M = f[0].length;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				target[i][j] = f[i][j] - x[i][j];
	}
	public static double[][] minus(double[][] f, double[][][] z, int thirdindex)
	{
		int N = f.length, M = f[0].length;
		double[][] a = new double [N][M];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				a[i][j] = f[i][j] - z[i][j][thirdindex];
		return a;
	}
	public static void log(double[][] f)
	{
		for (int i = 0; i < f.length; i++)
			for(int j = 0; j < f[0].length; j++)
				f[i][j] = Math.log(f[i][j]);
	}
	public static void log(double[] f)
	{
		for (int i = 0; i < f.length; i++)
				f[i] = Math.log(f[i]);
	}
	
	/**
	 * replaces the magnitude of each complex number by it's logarithm, scaling so that everything
	 * is positive.
	 * @param f
	 * @param min
	 */
	public static void logComplexMagRescale(double[][][] f, double fixmin)
	{
		double min = Math.max(FieldOps.magMin(f), fixmin);
		double logmag;
		double x, y;
		for (int i = 0; i < f.length; i++)
			for(int j = 0; j < f[0].length; j++)
			{
				f[i][j][0] *= 1/min;
				f[i][j][1] *= 1/min;
				if (Complex.mag(f[i][j]) < 1)
				{
					f[i][j][0] = 1;
					f[i][j][1] = 0;
				}
				logmag = Math.log(Complex.mag(f[i][j]));
				x = logmag*Math.cos(Complex.phase(f[i][j]));
				y = logmag*Math.sin(Complex.phase(f[i][j]));
				f[i][j][0] = x;
				f[i][j][1] = y;
			}
	}
	public static void negate(double[][] f)
	{
		for (int i = 0; i < f.length; i++)
			for(int j = 0; j < f[0].length; j++)
				f[i][j] = -f[i][j];
	}
	public static void negate(boolean[][] f)
	{
		for (int i = 0; i < f.length; i++)
			for(int j = 0; j < f[0].length; j++)
				f[i][j] = !f[i][j];
	}
	public static double[][] getLog(double[][] f)
	{
		double[][] l = new double [f.length][f[0].length];
		for (int i = 0; i < f.length; i++)
			for(int j = 0; j < f[0].length; j++)
				l[i][j] = Math.log(f[i][j]);
		return l;
	}
	public static double magMax(double[][][] f)
	{
		double max = -Double.MAX_VALUE;
		for (int i = 0; i < f.length; i++)
			for(int j = 0; j < f[0].length; j++)
				max = Math.max(max, Complex.mag(f[i][j]));
		return max;
	}
	public static double magMax(double[][] x, double[][] y)
	{
		double max = -Double.MAX_VALUE;
		for (int i = 0; i < x.length; i++)
			for(int j = 0; j < x[0].length; j++)
				max = Math.max(max, Distance.distance(x[i][j], y[i][j]));
		return max;
	}
	public static double magMin(double[][][] f)
	{
		double min = Double.MAX_VALUE;
		for (int i = 0; i < f.length; i++)
			for(int j = 0; j < f[0].length; j++)
				min = Math.min(min, Complex.mag(f[i][j]));
		return min;
	}
	
	public static double atan(double x, double y)
	{
		if (x == 0) return y > 0 ? Math.PI/2 : 3*Math.PI/2;
		if (x < 0) return Math.atan(y/x) + Math.PI;
		if (y < 0) return Math.atan(y/x) + 2*Math.PI;
		return Math.atan(y/x);
	}
	/**
	 * This produces angles between -pi and pi, again assuming zero = positive x axis.
	 * @param x
	 * @param y
	 * @return
	 */
	public static double atanCenteredZero(double x, double y)
	{
		if (x == 0) return y > 0 ? Math.PI/2 : -Math.PI/2;
		if (x > 0) return Math.atan(y/x);
		if (y > 0) return Math.atan(y/x) + Math.PI;
		if (y < 0) return Math.atan(y/x) - Math.PI;
		if (y == 0) return Math.PI;
		return 0;
	}
	
	//the weighting function is f^2.
	public static double sigmaX(double[][] f, int xi, int xf, int yi, int yf)
	{
		double var = 0;
		double sum = 0;
		double mean = 0;
		for (int i = xi; i < xf; i++)
			for (int j = yi; j < yf; j++){
				var += i*i*f[i][j]*f[i][j];
				mean += i*f[i][j]*f[i][j];
				sum += f[i][j]*f[i][j];
			}
		var /= sum;
		mean /= sum;
		var -= mean*mean;
		return Math.sqrt(var);
	}
	public static double sigmaY(double[][] f, int xi, int xf, int yi, int yf)
	{
		double var = 0;
		double sum = 0;
		double mean = 0;
		for (int i = xi; i < xf; i++)
			for (int j = yi; j < yf; j++){
				var += j*j*f[i][j]*f[i][j];
				mean += j*f[i][j]*f[i][j];
				sum += f[i][j]*f[i][j];
			}
		var /= sum;
		mean /= sum;
		var -= mean*mean;
		return Math.sqrt(var);
	}
	public static double[] sigmaR(double[][] f)
	{
		return new double[] {sigmaX(f, 0, f.length, 0, f[0].length), sigmaY(f, 0, f.length, 0, f[0].length)};
	}
	
	public static double centX(double[][] f, int xi, int xf, int yi, int yf)
	{
		double c = 0;
		double sum  = 0;
		for (int i = xi; i < xf; i++)
			for (int j = yi; j < yf; j++){
				c += i*f[i][j];
				sum += f[i][j];
			}
		return c/sum;
	}
	public static double[] centroid(double[][] f, boolean[][] include)
	{
		double cx = 0, cy = 0;
		double sum  = 0;
		for (int i = 0; i < f.length; i++)
			for (int j = 0; j < f[i].length; j++)
				if (include[i][j]){
					cx += i*f[i][j];
					cy += j*f[i][j];
					sum += f[i][j];
				}
			
		return new double[] {cx/sum, cy/sum};
	}
	public static double centXSq(double[][] f, int xi, int xf, int yi, int yf)
	{
		double c = 0;
		double sum  = 0;
		for (int i = xi; i < xf; i++)
			for (int j = yi; j < yf; j++){
				c += i*f[i][j]*f[i][j];
				sum += f[i][j]*f[i][j];
			}
		return c/sum;
	}
	public static double centXSq(double[][] f)
	{
		return centXSq(f, 0, f.length, 0, f[0].length);
	}
	public static double centY(double[][] f, int xi, int xf, int yi, int yf)
	{
		double c = 0;
		double sum  = 0;
		for (int i = xi; i < xf; i++)
			for (int j = yi; j < yf; j++){
				c += j*f[i][j];
				sum += f[i][j];
			}
		return c/sum;
	}
	public static double centYSq(double[][] f, int xi, int xf, int yi, int yf)
	{
		double c = 0;
		double sum  = 0;
		for (int i = xi; i < xf; i++)
			for (int j = yi; j < yf; j++){
				c += j*f[i][j]*f[i][j];
				sum += f[i][j]*f[i][j];
			}
		return c/sum;
	}
	public static double centYSq(double[][] f)
	{
		return centYSq(f, 0, f.length, 0, f[0].length);
	}
	
	public static double mean(double[][] f)
	{
		double sum  = 0;
		for (int i = 0; i < f.length; i++)
			for (int j = 0; j < f[0].length; j++){
				sum += f[i][j];
			}
		return sum/(f.length*f[0].length);
	}
	public static double mean(double[][] f, boolean[][] selectedPoints)
	{
		double sum  = 0;
		int n = 0;
		for (int i = 0; i < f.length; i++)
			for (int j = 0; j < f[0].length; j++)
				if (selectedPoints[i][j]){
					sum += f[i][j];
					n++;
				}
		if (n == 0) return FieldOps.mean(f);
		else return sum/n;
	}
	public static double mean(double[][] f, double[][] weight)
	{
		double sum  = 0;
		double n = 0;
		for (int i = 0; i < f.length; i++)
			for (int j = 0; j < f[0].length; j++){
					sum += f[i][j]*weight[i][j];
					n += weight[i][j];
				}
		if (n == 0) return FieldOps.mean(f);
		else return sum/n;
	}
	public static double sum(double[][] f)
	{
		double sum  = 0;
		for (int i = 0; i < f.length; i++)
			for (int j = 0; j < f[0].length; j++)
					sum += f[i][j];
		
		return sum;
	}
	public static double mean(double[][] f, int[][] selectedPoints, int value)
	{
		double sum  = 0;
		int n = 0;
		for (int i = 0; i < f.length; i++)
			for (int j = 0; j < f[0].length; j++)
				if (selectedPoints[i][j] == value){
					sum += f[i][j];
					n++;
				}
		if (n == 0) return 0;
		else return sum/n;
	}
	public static double mean(double[][] f, int xi, int xf, int yi, int yf)
	{
		double sum  = 0;
		int imin, imax, jmin, jmax;
		imin = Math.max(xi, 0);
		imax = Math.min(xf, f.length);
		jmin = Math.max(yi, 0);
		jmax = Math.min(yf, f[0].length);
		int n = 0;
		for (int i = imin; i < imax; i++)
			for (int j = jmin; j < jmax; j++){
				sum += f[i][j];
				n++;
			}
		return sum/n;
	}
	public static double mean(double[][][] f, int xi, int xf, int yi, int yf, int thirdindex)
	{
		double sum  = 0;
		for (int i = xi; i < xf; i++)
			for (int j = yi; j < yf; j++){
				sum += f[i][j][thirdindex];
			}
		return sum/((xf-xi)*(yf-yi));
	}
	public static double sigma(double[][] f, int xi, int xf, int yi, int yf, double mean)
	{
		double sum  = 0;
		for (int i = xi; i < xf; i++)
			for (int j = yi; j < yf; j++){
				sum += f[i][j]*f[i][j];
			}
		sum /= ((xf-xi)*(yf-yi));
		sum -= mean*mean;
		return Math.sqrt(sum);
	}
	public static double sigma(double[][] f, double mean)
	{
		double sum  = 0;
		for (int i = 0; i < f.length; i++)
			for (int j = 0; j < f[0].length; j++){
				sum += f[i][j]*f[i][j];
			}
		sum /= ((f.length*f[0].length));
		sum -= mean*mean;
		return Math.sqrt(sum);
	}
	public static double sigmaLocal(double[][] f, double mean, int io, int jo, int di, int dj)
	{
		double sum  = 0;
		int ninsum = 0;
		for (int i = io; i < io+di; i++)
			for (int j = jo; j < jo+dj; j++){
				if (i < f.length && j < f[0].length && i >= 0 && j >= 0)
					sum += f[i][j]*f[i][j];
					ninsum++;
			}
		sum /= di*dj;
		sum -= mean*mean;
		return Math.sqrt(sum);
	}
	public static double sigma(double[][] f)
	{
		return sigma(f, mean(f));
	}
	public static void shift(double[][][] source, double[][][] target)
	{
		int N = target.length, M = target[0].length;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
			{
				target[i][j][0] = source[(i + N/2) % N][(j + M/2) % M][0];
				target[i][j][1] = source[(i + N/2) % N][(j + M/2) % M][1];
			}
	}
	public static void shift(boolean[][] source, boolean[][] target)
	{
		int N = target.length;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				target[i][j] = source[(i + N/2) % N][(j + N/2) % N];
			}
	}
	public static void shift(boolean[][] source, boolean[][] target, int dx, int dy)
	{
		int N = target.length;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				target[i][j] = source[(i + dx + N) % N][(j + dy + N) % N];
			}
	}
	/**
	 * Periodic boundary conditions. It is assumed that the data is square and
	 * that abs(dx,dy) < N
	 * @param source
	 * @param target
	 * @param dx
	 * @param dy
	 */
	public static void shift(double[][][] source, double[][][] target, int dx, int dy)
	{
		int N = target.length;
		int M = target[0].length;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
			{
				target[i][j] = source[(i + dx + N) % N][(j + dy + M) % M];
			}
	}
	public static void shift(double[][] source, double[][] target, int dx, int dy)
	{
		int N = target.length;
		int M = target[0].length;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
			{
				target[i][j] = source[(i + dx + N) % N][(j + dy + M) % M];
			}
	}
	public static void shiftAvgExcluded(double[][] source, double[][] target, int dx, int dy)
	{
		int N = target.length;
		int M = target[0].length;
		double mean = FieldOps.mean(source);
		
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
			{
				int p = i-dx;
				int q = j-dy;
				if (withinBounds(p, q, source))
					target[i][j] = source[p][q];
				else target[i][j] = mean;
			}
	}
	public static void shift(double[][] source, double[][] target)
	{
		int N = target.length;
		int M = target[0].length;
		shift(source, target, N/2, M/2);
	}

	public static double max(double[][] f, int xi, int xf, int yi, int yf)
	{
		double max = -Double.MAX_VALUE;
		for (int i = xi; i < xf; i++)
			for (int j = yi; j < yf; j++){
				if (f[i][j] != Double.NaN) 
					max = Math.max(max, f[i][j]);
			}
		return max;
	}
	public static double max(double[][] f, boolean[][] include)
	{
		double max = -Double.MAX_VALUE;
		for (int i = 0; i < f.length; i++)
			for (int j = 0; j < f[i].length; j++){
				if (include[i][j]) 
					max = Math.max(max, f[i][j]);
			}
		return max;
	}
	/**
	 * This returns false if all surrounding pixels have the same value as this pixel, true otherwise.
	 * @param data
	 * @return
	 */
	public static boolean[][] isBoundary(boolean[][] data){
		int nx = data.length, ny = data[0].length;
		boolean[][] ans = new boolean [nx][ny];
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				ans[i][j] = ans[i][j] || (withinBounds(i-1, j-1, data) && (data[i-1][j-1] != data[i][j]));
				ans[i][j] = ans[i][j] || (withinBounds(i+1, j-1, data) && (data[i+1][j-1] != data[i][j]));
				ans[i][j] = ans[i][j] || (withinBounds(i-1, j+1, data) && (data[i-1][j+1] != data[i][j]));
				ans[i][j] = ans[i][j] || (withinBounds(i+1, j+1, data) && (data[i+1][j+1] != data[i][j]));
				ans[i][j] = ans[i][j] || (withinBounds(i, j+1, data) && (data[i][j+1] != data[i][j]));
				ans[i][j] = ans[i][j] || (withinBounds(i, j-1, data) && (data[i][j-1] != data[i][j]));
				ans[i][j] = ans[i][j] || (withinBounds(i-1, j, data) && (data[i-1][j] != data[i][j]));
				ans[i][j] = ans[i][j] || (withinBounds(i+1, j, data) && (data[i+1][j] != data[i][j]));
			}
		return ans;
	}
	
	/**
	 * Code: 0 = minimum, 1 = maximum, 2 = either.
	 * @param data
	 * @param direction
	 * @param code
	 * @return
	 */
	public static boolean[][] isDirectionalExtremum(double[][] data, double[] direction, int code){
		double[][] directionalDerivative = directionalDerivative2D(data, direction);
//		LayerViewer.show(Layer.getFreeLayer(directionalDerivative), 1000, true);
		double[][] directional2Derivative = directionalDerivative2D(directionalDerivative, direction);
//		LayerViewer.show(Layer.getFreeLayer(directional2Derivative), 1000, true);
		boolean[][] tempb = FieldOps.isGreaterThan(directionalDerivative, 0);
//		LayerViewer.show(Layer.getFreeLayer(ArrayOps.toDouble(tempb)), 1000, true);
		boolean[][] zerosD = FieldOps.isBoundary(tempb);
//		LayerViewer.show(Layer.getFreeLayer(ArrayOps.toDouble(zerosD)), 1000, true);
		if (code == 2) return zerosD;
		boolean[][] ans = new boolean[data.length][data[0].length];
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[i].length; j++)
			{
				if (code == 0)
					ans[i][j] = zerosD[i][j] && directional2Derivative[i][j] > 0;
				else if (code == 1)
					ans[i][j] = zerosD[i][j] && directional2Derivative[i][j] < 0;
			}
		return ans;
	}
	public static boolean[][] isDirectionalExtremum(double[][] data, double[] direction, int code, double twoDcutoff){
		double[][] directionalDerivative = directionalDerivative2D(data, direction);
		LayerViewer.show(Layer.getFreeLayer(directionalDerivative), 1000, true);
		double[][] directional2Derivative = directionalDerivative2D(directionalDerivative, direction);
		LayerViewer.show(Layer.getFreeLayer(directional2Derivative), 1000, true);
		boolean[][] tempb = FieldOps.isGreaterThan(directionalDerivative, 0);
		LayerViewer.show(Layer.getFreeLayer(ArrayOps.toDouble(tempb)), 1000, true);
		boolean[][] zerosD = FieldOps.isBoundary(tempb);
		LayerViewer.show(Layer.getFreeLayer(ArrayOps.toDouble(zerosD)), 1000, true);
		if (code == 2) return zerosD;
		boolean[][] ans = new boolean[data.length][data[0].length];
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[i].length; j++)
			{
				if (code == 0)
					ans[i][j] = zerosD[i][j] && directional2Derivative[i][j] > twoDcutoff;
				else if (code == 1)
					ans[i][j] = zerosD[i][j] && directional2Derivative[i][j] < twoDcutoff;
			}
		return ans;
	}
	public static double[] getSpectrumWithBools(double[][] data, boolean[][][] filter)
	{
//		String outdir = dir + ending + "_" + cutname + "_cut\\";// + (normalizeEachLine ? "norm_layer" : "") + "\\";
//		File outd = new File (outdir);
//		if (!outd.exists()) outd.mkdir();
		
					
		int npts = filter.length;
//		double[][] onData = stripCutEachEnergy(t, width, unitBragg.getA());
		double[] result = new double [npts];
		for (int i = 0; i < npts; i++)
			result[i] = FieldOps.mean(data, filter[i]);
		return result;
//		double[][] resultTranspose = FieldOps.transpose(result);
		//obtain E vs K. For K, we remember that the dataset is the symmetrized FFT so that it's real space is actually k-space.
		
//		PointSpectra on = new PointSpectra(resultTranspose, xAxis, new double[t.nlayers], t.v);
//		PointSpectra.writeBIN(on, outdir + "Intensity spectra I vs k spec.bin");
//		Layer onl = on.toLayer();
//		onl.flipY();
//		Layer.writeBIN(onl, outdir + "Intensity spectra I vs k layer.bin");
		
//		if (JOptionPane.showConfirmDialog(null, "Write the masks used?") == JOptionPane.YES_OPTION)
//			Topomap.writeBIN(new Topomap(ArrayOps.toDouble(filter), xAxis, t.x, t.y, null), outdir + "masks.bin");
//		//
	}
	
	public static boolean[][][] isDirectionalExtremum(double[][][] data, double[] direction, int code){
		boolean[][][] ans = new boolean[data.length][][];
		for (int i = 0; i < data.length; i++)
			ans[i] = isDirectionalExtremum(data[i], direction, code);
		return ans;
	}
	public static boolean[][][] isDirectionalExtremum(double[][][] data, double[] direction, int code, double twoDcutoff){
		boolean[][][] ans = new boolean[data.length][][];
		for (int i = 0; i < data.length; i++)
			ans[i] = isDirectionalExtremum(data[i], direction, code, twoDcutoff);
		return ans;
	}

	public static int max(int[][] f, int xi, int xf, int yi, int yf)
	{
		int max = Integer.MIN_VALUE;
		for (int i = xi; i < xf; i++)
			for (int j = yi; j < yf; j++){
				max = Math.max(max, f[i][j]);
			}
		return max;
	}
	public static double max(double[][] f)
	{
		return max(f, 0, f.length, 0, f[0].length);
	}
	public static int max(int[][] f)
	{
		return max(f, 0, f.length, 0, f[0].length);
	}
	public static double min(double[][] f, int xi, int xf, int yi, int yf)
	{
		double min = Double.MAX_VALUE;
		for (int i = xi; i < xf; i++)
			for (int j = yi; j < yf; j++){
				if (f[i][j] != Double.NaN) 
					min = Math.min(min, f[i][j]);
			}
		return min;
	}
	public static double min(double[][][] f)
	{
		double min = Double.MAX_VALUE;
		for (int i = 0; i < f.length; i++)
			for (int j = 0; j < f[i].length; j++)
				for (int k = 0; k < f[i][j].length; k++)
					if (f[i][j][k] != Double.NaN) 
						min = Math.min(min, f[i][j][k]);
		return min;
	}
	public static double min(double[][] f, boolean[][] include)
	{
		double min = Double.MAX_VALUE;
		for (int i = 0; i < f.length; i++)
			for (int j = 0; j < f[i].length; j++){
				if (include[i][j]) 
					min = Math.min(min, f[i][j]);
			}
		return min;
	}

	public static int min(int[][] f, int xi, int xf, int yi, int yf)
	{
		int min = Integer.MAX_VALUE;
		for (int i = xi; i < xf; i++)
			for (int j = yi; j < yf; j++){
				min = Math.min(min, f[i][j]);
			}
		return min;
	}
	public static double min(double[][] f)
	{
		return min(f, 0, f.length, 0, f[0].length);
	}
	public static int min(int[][] f)
	{
		return min(f, 0, f.length, 0, f[0].length);
	}
	public static int[] minR(double[][] f, int xi, int xf, int yi, int yf)
	{
		int[] minR = new int[2];
		double min = Double.MAX_VALUE;
		double oldmin = min;
		for (int i = xi; i < xf; i++)
			for (int j = yi; j < yf; j++){
				min = Math.min(min, f[i][j]);
				if (min != oldmin){
					minR[0] = i; minR[1] = j;
					oldmin = min;
				}
			}
		return minR;
	}
	public static int[] minR(double[][] f)
	{
		int[] minR = new int[2];
		double min = Double.MAX_VALUE;
		double oldmin = min;
		for (int i = 0; i < f.length; i++)
			for (int j = 0; j < f[i].length; j++){
				min = Math.min(min, f[i][j]);
				if (min != oldmin){
					minR[0] = i; minR[1] = j;
					oldmin = min;
				}
			}
		return minR;
	}
	public static int[] maxR(double[][] f, int xi, int xf, int yi, int yf)
	{
		int[] maxR = new int[2];
		double max = -Double.MAX_VALUE;
		double oldmax = max;
		for (int i = xi; i < xf; i++)
			for (int j = yi; j < yf; j++){
				max = Math.max(max, f[i][j]);
				if (max != oldmax){
					maxR[0] = i; maxR[1] = j;
					oldmax = max;
				}
			}
		return maxR;
	}
	
	//This returns a field of twice the size of the input field (a square), with values linearly interpolated between
	//the values of the input field.
	public static double[][] expandNoSmooth(double[][] in, int factor)
	{
		int N = in.length*factor;
		double[][] out = new double[N][N];
		for (int i =0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				out[i][j] = in[i/factor][j/factor];
			}
			
		return out;
	}
	public static double[][] expand(double[][] in)
	{
		int N = in.length*2;
		double[][] out = new double[N][N];
		for (int i = 0; i < N/2; i++)
			for (int j = 0; j < N/2; j++)
			{
				out[2*i][2*j] = in[i][j];
			}
		//interpolate the 4 edges
		for (int i = 0; i < N/2-1; i++)
			out[2*i + 1][0] = (in[i][0] + in[i+1][0])/2;
		for (int i = 0; i < N/2-1; i++)
			out[0][2*i + 1] = (in[0][i] + in[0][i+1])/2;
		for (int i = 0; i < N/2-1; i++)
			out[2*i + 1][N-1] = (in[i][(N/2)-1] + in[i+1][(N/2)-1])/2;
		for (int i = 0; i < N/2-1; i++)
			out[N-1][2*i + 1] = (in[(N/2)-1][i] + in[(N/2)-1][i+1])/2;
		for (int i = 0; i < N/2-1; i++)
			out[2*i][N-1] = in[i][(N/2)-1];
		for (int i = 0; i < N/2-1; i++)
			out[N-1][2*i] = in[(N/2)-1][i];
		
		//the corner pieces must be...
		out[N-2][N-1] = in[(N/2)-1][(N/2)-1];
		out[N-1][N-2] = in[(N/2)-1][(N/2)-1];
		out[N-1][N-1] = in[(N/2)-1][(N/2)-1];
		
		//if both are even, it came from in. If one is even, average the two
		//nearest neighbors. If both are odd, average the four corner neighbors.
		int ini, inj;
		for (int i = 1; i < N-1; i++)
			for (int j = 1; j < N-1; j++)
			{
				ini = i/2; inj = j/2;
				if (i%2 != 0 && j%2 != 0)
					out[i][j] = (in[ini][inj] + in[ini+1][inj] + in[ini][inj+1] + in[ini+1][inj+1])/4;
				else if (i%2 == 0 && j%2 != 0)
					out[i][j] = (in[ini][inj] + in[ini][inj+1])/2;
				else if (j%2 == 0 && i%2 != 0)
					out[i][j] = (in[ini][inj] + in[ini+1][inj])/2;
			}
			
		return out;
	}
	
	/**
	 * This attempts to expand something by adding a margin to the Fourier Transform, then inverse Fourier-transforming that.
	 * @param in
	 * @param n
	 * @return
	 */
	public static double[][] expandFourier(double[][] in, int nx, int ny)
	{
		int lx = in.length, ly = in[0].length;
		int marginX = nx-lx, marginY = ny-ly;
		
		double[][][] fftz2 = FFTOps.obtainFFT(in).fHat2();
		double magmin = magMin(fftz2);
		double[][][] newFFTZ2 = new double [nx][ny][2];
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny/2; j++)
			{
				double theta = 2*Math.PI*Math.random();
				newFFTZ2[i][j][0] = magmin*Math.cos(theta);
				newFFTZ2[i][j][1] = magmin*Math.sin(theta);
				newFFTZ2[nx-i-1][ny-j-1][0] = magmin*Math.cos(theta);
				newFFTZ2[nx-i-1][ny-j-1][1] = -magmin*Math.sin(theta);
			}
		for (int i = marginX/2; i < (marginX/2)+lx; i++)
			for (int j = marginY/2; j < (marginY/2)+ly; j++)
			{
				newFFTZ2[i][j][0] = fftz2[i-marginX/2][j-marginY/2][0];
				newFFTZ2[i][j][1] = fftz2[i-marginX/2][j-marginY/2][1];
			}
		
		double[][] ans = new double [nx][ny];
		FFTOps.obtainIFFTCent(newFFTZ2, ans);
		return ans;
	}
	
	//expands the portion of in between imin, imax, jmin, jmax
	public static double[][] expand(double[][] in, int imin, int imax, int jmin, int jmax)
	{
		int N = (imax - imin)*2;
		int M = (jmax - jmin)*2;
		double[][] out = new double[N][M];
		for (int i = 0; i < N/2; i++)
			for (int j = 0; j < M/2; j++)
			{
				out[2*i][2*j] = in[imin+i][jmin+j];
			}
		//interpolate the 4 edges
		for (int i = 0; i < N/2-1; i++)
			out[2*i + 1][0] = (in[imin+i][jmin+0] + in[imin+i+1][jmin+0])/2;
		for (int i = 0; i < N/2-1; i++)
			out[0][2*i + 1] = (in[imin+0][jmin+i] + in[imin+0][jmin+i+1])/2;
		for (int i = 0; i < N/2-1; i++)
			out[2*i + 1][N-1] = (in[imin+i][jmin+(N/2)-1] + in[imin+i+1][jmin+(N/2)-1])/2;
		for (int i = 0; i < N/2-1; i++)
			out[N-1][2*i + 1] = (in[imin+(N/2)-1][jmin+i] + in[imin+(N/2)-1][jmin+i+1])/2;
		for (int i = 0; i < N/2-1; i++)
			out[2*i][N-1] = in[imin+i][jmin+(N/2)-1];
		for (int i = 0; i < N/2-1; i++)
			out[N-1][2*i] = in[imin+(N/2)-1][jmin+i];
		
		//the corner pieces must be...
		out[N-2][N-1] = in[imin+(N/2)-1][jmin+(N/2)-1];
		out[N-1][N-2] = in[imin+(N/2)-1][jmin+(N/2)-1];
		out[N-1][N-1] = in[imin+(N/2)-1][jmin+(N/2)-1];
		
		//if both are even, it came from in. If one is even, average the two
		//nearest neighbors. If both are odd, average the four corner neighbors.
		int ini, inj;
		for (int i = 1; i < N-1; i++)
			for (int j = 1; j < N-1; j++)
			{
				ini = i/2; inj = j/2;
				if (i%2 != 0 && j%2 != 0)
					out[i][j] = (in[imin+ini][jmin+inj] + in[imin+ini+1][jmin+inj] + in[imin+ini][jmin+inj+1] + in[imin+ini+1][jmin+inj+1])/4;
				else if (i%2 == 0 && j%2 != 0)
					out[i][j] = (in[imin+ini][jmin+inj] + in[imin+ini][jmin+inj+1])/2;
				else if (j%2 == 0 && i%2 != 0)
					out[i][j] = (in[imin+ini][jmin+inj] + in[imin+ini+1][jmin+inj])/2;
			}
			
		return out;
	}
	public static double[][] expand4(double[][] in)
	{
		return expand(expand(in));
	}
	public static double[][] expand4(double[][] in, int imin, int imax, int jmin, int jmax)
	{
		return expand(expand(in, imin, imax, jmin, jmax));
	}
	
	public static double[][] zeropad(double[][] in)
	{
		double[][] ans = new double [in.length*2][in[0].length*2];
		for (int i = 0; i < in.length; i++)
			for (int j = 0; j < in[0].length; j++)
				ans[i][j] = in[i][j];
		return ans;
	}
	
	public static double[][][] expand(double[][][] in)
	{
		int N = in.length*2;
		double[][][] out = new double[N][N][2];
		for (int k = 0; k < 2; k++){
			for (int i = 0; i < N/2; i++)
				for (int j = 0; j < N/2; j++)
				{
					out[2*i][2*j] = in[i][j];
				}
			//interpolate the 4 edges
			for (int i = 0; i < N/2-1; i++)
				out[2*i + 1][0][k] = (in[i][0][k] + in[i+1][0][k])/2;
			for (int i = 0; i < N/2-1; i++)
				out[0][2*i + 1][k] = (in[0][i][k] + in[0][i+1][k])/2;
			for (int i = 0; i < N/2-1; i++)
				out[2*i + 1][N-1][k] = (in[i][(N/2)-1][k] + in[i+1][(N/2)-1][k])/2;
			for (int i = 0; i < N/2-1; i++)
				out[N-1][2*i + 1][k] = (in[(N/2)-1][i][k] + in[(N/2)-1][i+1][k])/2;
			for (int i = 0; i < N/2-1; i++)
				out[2*i][N-1][k] = in[i][(N/2)-1][k];
			for (int i = 0; i < N/2-1; i++)
				out[N-1][2*i][k] = in[(N/2)-1][i][k];
			
			//the corner pieces must be...
			out[N-2][N-1][k] = in[(N/2)-1][(N/2)-1][k];
			out[N-1][N-2][k] = in[(N/2)-1][(N/2)-1][k];
			out[N-1][N-1][k] = in[(N/2)-1][(N/2)-1][k];
			
			//if both are even, it came from in. If one is even, average the two
			//nearest neighbors. If both are odd, average the four corner neighbors.
			int ini, inj;
			for (int i = 1; i < N-1; i++)
				for (int j = 1; j < N-1; j++)
				{
					ini = i/2; inj = j/2;
					if (i%2 != 0 && j%2 != 0)
						out[i][j][k] = (in[ini][inj][k] + in[ini+1][inj][k] + in[ini][inj+1][k] + in[ini+1][inj+1][k])/4;
					else if (i%2 == 0 && j%2 != 0)
						out[i][j][k] = (in[ini][inj][k] + in[ini][inj+1][k])/2;
					else if (j%2 == 0 && i%2 != 0)
						out[i][j][k] = (in[ini][inj][k] + in[ini+1][inj][k])/2;
				}
			}
		return out;
	}
	
	/**
	 * This subtracts from data a linear function of the form c1*x + c2*y, 
	 * where c1 and c2 are determined (in the y-case) c2 = (average of top half - average of bottom half)/distance between their centers.
	 * @param data
	 */
	public static void subtractSlopes1(double[][] data)
	{
		int N = data.length; //square table
		double qss = 0, qsm = 0, qms = 0, qmm = 0;
		for (int i = 0; i < N/2; i++)
			for (int j = 0; j < N/2; j++)
				qss += data[i][j];
		for (int i = N/2; i < N; i++)
			for (int j = 0; j < N/2; j++)
				qms += data[i][j];
		for (int i = 0; i < N/2; i++)
			for (int j = N/2; j < N; j++)
				qsm += data[i][j];
		for (int i = N/2; i < N; i++)
			for (int j = N/2; j < N; j++)
				qmm += data[i][j];
		
		qss /= N*N/4;
		qms /= N*N/4;
		qsm /= N*N/4;
		qmm /= N*N/4;
		double hsx, hsy, hmx, hmy;
		hsx = (qss+qsm)/2;
		hsy = (qss+qms)/2;
		hmx = (qmm+qms)/2;
		hmy = (qmm+qsm)/2;
		double sx, sy;
		sx = 2*(hmx-hsx)/N;
		sy = 2*(hmy-hsy)/N;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				data[i][j] -= sx*i + sy*j;
		
	}
	public static void subtractAvg(double[][] data)
	{
		int N = data.length; 
		int M = data[0].length;//square table
		double mean = FieldOps.mean(data, 0, N, 0, M);
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				data[i][j] -= mean;
	}
	public static void subtractAvgFactor(double[][] data, double factor)
	{
		double mean = FieldOps.mean(data);
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
				data[i][j] -= mean*factor;
	}
	public static double[][] subset(double[][] data, int i0, int i1, int j0, int j1)
	{
		double[][] ans = new double [i1-i0][j1-j0];
		for (int i = i0; i < i1; i++)
			for (int j = j0; j < j1; j++)
				ans[i-i0][j-j0] = data[i][j];
		return ans;
	}
	/**
	 * This attempts to return the "constant height dI/dV." Assuming that for the entire range of Z,
	 * dI/dV(x,y) is proportional to exp(-kZ(x,y)), we multiply the dI/dV by exp(kZ(x,y)) to extract what the 
	 * dI/dV "would have been" at constant height. See e.g. PRL 81, 5616 (1998).
	 * 
	 * To keep the order of magnitude the same, Z has its minimum value subtracted before exponentiation.
	 * @param didv
	 * @param z
	 * @param k
	 * @return
	 */
	public static double[][] getZAdjusteddIdV(double[][] didv, double[][] z, double k, boolean renormalize)
	{
		int nx = z.length, ny = z[0].length;
		double[][] ans = new double [nx][ny];
		double[][] z0 = null;
		if (renormalize) z0 = FieldOps.minus(z, FieldOps.min(z));
		else z0 = z;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				ans[i][j] = didv[i][j]*Math.exp(k*z0[i][j]);
		return ans;
	}
	
	/**
	 * This returns the magnitude of the smoothed function exp(iq dot r)*data with smoothing length l.
	 * 
	 * This is basically the same as the u-field calculation except we take the magnitude instead of the phase.
	 * @param data
	 * @param q
	 * @param lPixels
	 * @return
	 */
	public static double[][] getSignalAmplitude(double[][] data, double[] q, double l)
	{
		double mean = FieldOps.mean(data);
		FieldOps.minusEquals(data, mean);
		int nx = data.length, ny = data[0].length;
		int N = nx;
		double[] braggTrue = new double[2];
		for (int i = 0; i < 2; i++)
				braggTrue[i] = q[i]*2*Math.PI/N;

		double[][][] cossin = new double [2][nx][ny];
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				cossin[0][i][j] = Math.cos(dot(braggTrue, i, j))*data[i][j];
				cossin[1][i][j] = -Math.sin(dot(braggTrue, i, j))*data[i][j];
			}
		
		double[][] cosSmooth = FieldOps.gaussSmooth(cossin[0], l);
		double[][] sinSmooth = FieldOps.gaussSmooth(cossin[1], l);
		FieldOps.plusEquals(data, mean);
		return FieldOps.magnitude(cosSmooth, sinSmooth);
		
		
	}
	public static double[][][] getSignalZ_2NM(double[][] data, double[] q, double l)
	{
		double mean = FieldOps.mean(data);
		FieldOps.minusEquals(data, mean);
		int nx = data.length, ny = data[0].length;
		int N = nx;
		double[] braggTrue = new double[2];
		for (int i = 0; i < 2; i++)
				braggTrue[i] = q[i]*2*Math.PI/N;

		double[][][] cossin = new double [2][nx][ny];
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				cossin[0][i][j] = Math.cos(dot(braggTrue, i, j))*data[i][j];
				cossin[1][i][j] = -Math.sin(dot(braggTrue, i, j))*data[i][j];
			}
		
		double[][] cosSmooth = FieldOps.gaussSmooth(cossin[0], l);
		double[][] sinSmooth = FieldOps.gaussSmooth(cossin[1], l);
		FieldOps.plusEquals(data, mean);
		return new double[][][] {cosSmooth, sinSmooth};
		
	}
	public static double dot(double[] a, int b0, int b1)
	{
		return a[0]*b0 + a[1]*b1;
	}

	public static void subtactLineAvg(double[][] data, boolean truei)
	{
		double mean = 0;
		int Ni = data.length;
		int Nj = data[0].length;
		if (truei)
			for (int i = 0; i < Ni; i++)
			{
				mean = ArrayOps.mean(data[i]);
				for (int j = 0; j < Nj; j++)
					data[i][j]-=mean;
			}
		else
			for (int i = 0; i < Nj; i++)
			{
				mean = 0;
				for (int j = 0; j < Ni; j++){
					mean += data[j][i];
				}
				mean /=Ni;
				for (int j = 0; j < Ni; j++){
					data[j][i] -= mean;
				}
			}
			
	}
	/**
	 * Subtracts the best fit line of each line in the square matrix data. If horizontal, the line is defined by keeping the 2nd index constant
	 * and sweeping the first; if !horizontal by keeping the 1st index constant and sweeping the 2nd.
	 * @param data
	 * @param truei
	 */
	public static void subtractLineFit(double[][] data, boolean horizontal)
	{
		double mean = 0;
		int Ni = data.length;
		int Nj = data[0].length;
		double[] line;
		double[] x;
		double[] p;
		Regression reg;
		if (horizontal)
		{
			line = new double [Ni];
			x = new double [Ni];
			for (int i = 0; i < Ni; i++)
				x[i] = i;
		}
		else
		{
			line = new double [Nj];
			x = new double [Nj];
			for (int j = 0; j < Nj; j++)
				x[j] = j;
		}
		
		
		
		if (horizontal)
			for (int j = 0; j < Nj; j++)
			{
				for (int i = 0; i < Ni; i++)
					line[i] = data[i][j];
				reg = new Regression(x, line);
				reg.linear();
				p = reg.getBestEstimates();
				for (int i = 0; i < Ni; i++)
					data[i][j] -= p[0] + p[1]*x[i];
			}
		else
			for (int i = 0; i < Ni; i++)
			{
				for (int j = 0; j < Nj; j++)
					line[j] = data[i][j];
				reg = new Regression(x, line);
				reg.linear();
				p = reg.getBestEstimates();
				for (int j = 0; j < Nj; j++)
					data[i][j] -= p[0] + p[1]*x[j];
			}
			
	}
	public static double subtractAvg(double[][] data, int x, int y, int di, int dj)
	{
		double mean = FieldOps.mean(data, x, x+di, y, y+dj);
		int imin, imax, jmin, jmax;
		imin = Math.max(x, 0);
		imax = Math.min(x+di, data.length);
		jmin = Math.max(y, 0);
		jmax = Math.min(y+dj, data[0].length);
		for (int i = imin; i < imax; i++)
			for (int j = jmin; j < jmax; j++)
				data[i][j] -= mean;
		return mean;
	}
	public static void addValue(double[][] data, int x, int y, int di, int dj, double value)
	{
		int imin, imax, jmin, jmax;
		imin = Math.max(x, 0);
		imax = Math.min(x+di, data.length);
		jmin = Math.max(y, 0);
		jmax = Math.min(y+dj, data[0].length);
		for (int i = imin; i < imax; i++)
			for (int j = jmin; j < jmax; j++)
				data[i][j] += value;
	}
	public static void subtractFracAvg(double[][] data, double frac)
	{
		int N = data.length; //square table
		double mean = FieldOps.mean(data, 0, N, 0, N);
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				data[i][j] -= frac*mean;
	}
	//This subtracts the average of each component of data, which is assuned to be 
	//a 2-component vector array of size [N][N][2].
	public static void subtractAvg(double[][][] data)
	{
		int N = data.length; //square table
		double[] mean = new double [2];
		for (int i = 0; i < 2; i++) mean[i] = FieldOps.mean(data, 0, N, 0, N, i);
		
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				for (int k = 0; k < 2; k++)
				data[i][j][k] -= mean[k];
	}
	//This returns an array indicating the number of times the 2pi line has been crossed in a 
	//discontinous "phase field." The output is an integer at each point. The addition 2pi*the integer
	//to the original field should render this field continuous.
	//The phase field is restricted to values between 0 and 2pi.
	public static int[][] phaseSteps(double[][] phase, int sx, int sy)
	{
		double p = Math.PI*2;
		int N = phase.length; //square
		int[][] n = new int[N][N];
		boolean[][] visited = new boolean [N][N];
		n[sx][sy] = 0;
		visited[sx][sy] = true;
		ArrayList<int[]> queue = new ArrayList<int[]>();
		queue.add(new int[] {sx, sy, sx+1, sy});
		queue.add(new int[] {sx, sy, sx-1, sy});
		queue.add(new int[] {sx, sy, sx, sy+1});
		queue.add(new int[] {sx, sy, sx, sy-1});
		int[] points;
		int counter = 0;
		while(queue.size() > 0)
		{
//			if (counter % 10000 == 0) System.out.println(counter + "\t" + queue.size());
			counter++;
			points = queue.remove(0);
			visited[points[2]][points[3]] = true;
			phaseStep(phase, points, p, n);
			if (points[2] > 0 && !visited[points[2]-1][points[3]]){
				queue.add(new int[] {points[2], points[3], points[2]-1, points[3]}); visited[points[2]-1][points[3]] = true;}
			if (points[3] > 0 && !visited[points[2]][points[3]-1]){
				queue.add(new int[] {points[2], points[3], points[2], points[3]-1}); visited[points[2]][points[3]-1] = true;}
			if (points[2] < N-1 && !visited[points[2]+1][points[3]]){
				queue.add(new int[] {points[2], points[3], points[2]+1, points[3]}); visited[points[2]+1][points[3]] = true;}
			if (points[3] < N-1 && !visited[points[2]][points[3]+1]){
				queue.add(new int[] {points[2], points[3], points[2], points[3]+1}); visited[points[2]][points[3]+1] = true;}
		}
		return n;
	}
	public static void putPhaseSteps(double[][] phase, int sx, int sy, int[][] n)
	{
		double p = Math.PI*2;
		int N = phase.length; //square
		boolean[][] visited = new boolean [N][N];
		n[sx][sy] = 0;
		visited[sx][sy] = true;
		ArrayList<int[]> queue = new ArrayList<int[]>();
		queue.add(new int[] {sx, sy, sx+1, sy});
		queue.add(new int[] {sx, sy, sx-1, sy});
		queue.add(new int[] {sx, sy, sx, sy+1});
		queue.add(new int[] {sx, sy, sx, sy-1});
		int[] points;
		int counter = 0;
		while(queue.size() > 0)
		{
//			if (counter % 10000 == 0) System.out.println(counter + "\t" + queue.size());
			counter++;
			points = queue.remove(0);
			visited[points[2]][points[3]] = true;
			phaseStep(phase, points, p, n);
			if (points[2] > 0 && !visited[points[2]-1][points[3]]){
				queue.add(new int[] {points[2], points[3], points[2]-1, points[3]}); visited[points[2]-1][points[3]] = true;}
			if (points[3] > 0 && !visited[points[2]][points[3]-1]){
				queue.add(new int[] {points[2], points[3], points[2], points[3]-1}); visited[points[2]][points[3]-1] = true;}
			if (points[2] < N-1 && !visited[points[2]+1][points[3]]){
				queue.add(new int[] {points[2], points[3], points[2]+1, points[3]}); visited[points[2]+1][points[3]] = true;}
			if (points[3] < N-1 && !visited[points[2]][points[3]+1]){
				queue.add(new int[] {points[2], points[3], points[2], points[3]+1}); visited[points[2]][points[3]+1] = true;}
		}
	}
	public static void putAddedPhaseSteps(double[][] phase, int[][] n, double p, double[][] target)
	{
		for (int i = 0; i < phase.length; i++)
			for (int j = 0; j < phase[0].length; j++)
				target[i][j] = phase[i][j] + n[i][j]*p;
	}
	
	//This method counts the number of disconnected "true" zones in a boolean map. The return value is an array of length
	//[3][n] where n is the number of such regions. the first entry is the x centroid, the second the y, the third the area of each region.
	
	public static double[][] isolatedTrueRegions(boolean[][] map)
	{
		double p = Math.PI*2;
		int N = map.length; //square
		boolean[][] visited = new boolean [N][N];
		ArrayList<int[]> queue = new ArrayList<int[]>();
		int[] points;
		int area = 0;
		double cx = 0, cy = 0;
		ArrayList<int[]> pointsinregion = new ArrayList<int[]>();
		ArrayList<Integer> areas = new ArrayList<Integer>();
		ArrayList<double[]> centroids = new ArrayList<double[]>();
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				if (map[i][j] && !visited[i][j])
				{
					area = 0; cx = 0; cy = 0;
					pointsinregion = new ArrayList<int[]>();
					queue.add(new int[] {i, j, i, j});
					//blot the region:
					while(queue.size() > 0)
					{
						points = queue.remove(0);
						pointsinregion.add(new int [] {points[2], points[3]});
						area++;
						visited[points[2]][points[3]] = true;
						if (points[2] > 0 && !visited[points[2]-1][points[3]] && map[points[2]-1][points[3]]){
							queue.add(new int[] {points[2], points[3], points[2]-1, points[3]}); visited[points[2]-1][points[3]] = true;}
						if (points[3] > 0 && !visited[points[2]][points[3]-1] && map[points[2]][points[3]-1]){
							queue.add(new int[] {points[2], points[3], points[2], points[3]-1}); visited[points[2]][points[3]-1] = true;}
						if (points[2] < N-1 && !visited[points[2]+1][points[3]] && map[points[2]+1][points[3]]){
							queue.add(new int[] {points[2], points[3], points[2]+1, points[3]}); visited[points[2]+1][points[3]] = true;}
						if (points[3] < N-1 && !visited[points[2]][points[3]+1] && map[points[2]][points[3]+1]){
							queue.add(new int[] {points[2], points[3], points[2], points[3]+1}); visited[points[2]][points[3]+1] = true;}
					}
					//record its data:
					for (int k = 0; k < pointsinregion.size(); k++)
					{
						cx += pointsinregion.get(k)[0];
						cy += pointsinregion.get(k)[1];
					}
					cx /= pointsinregion.size(); cy /= pointsinregion.size();
					centroids.add(new double[] {cx, cy});
					areas.add(area);
				}
			}
		//finally tally up the results:
		double[][] answer = new double [3][areas.size()];
		for (int i = 0; i < answer[0].length; i++)
		{
			answer[0][i] = centroids.get(i)[0]; 
			answer[1][i] = centroids.get(i)[1];
			answer[2][i] = areas.get(i);
		}
		
		return answer;
	}
	
	//returns the value of n at i2, j2, given its value at i, j.
	private static void phaseStep(double[][] phase, int[] points, double p, int[][] n)
	{
		//points is {i, j, i2, j2}
		//is it slowly varying?
		if (Math.abs(phase[points[0]][points[1]] - phase[points[2]][points[3]]) < p/4)
			n[points[2]][points[3]] = n[points[0]][points[1]];
		else if (phase[points[0]][points[1]] > phase[points[2]][points[3]])
			n[points[2]][points[3]] = n[points[0]][points[1]] + 1;
		else n[points[2]][points[3]] = n[points[0]][points[1]] - 1;
	}
	
	public static void translateOnePixel(double[][] source, double[][] target, int i, int j, double dx, double dy)
	{
		//first assume that the target position is completely within the range.
		int N = source.length;
		double v = source[i][j];
		int mpx = dx > 0 ? (int)(dx + 0.5) : (int)(dx - 0.5);
		int mpy = dy > 0 ? (int)(dy + 0.5) : (int)(dy - 0.5);
		double dxr = Math.abs(dx - mpx), dyr = Math.abs(dy - mpy);
		double sdx = 1-dxr, sdy = 1-dyr; 
		double small = dxr*dyr;
		double large = sdx*sdy;
		double xcut = dxr*sdy;
		double ycut = dyr*sdx;
//		if (dx - (int)dx < 0.01)
//		{	System.out.println(small + "\t" + large + "\t" + xcut + "\t" + ycut + "\t" + dxr + "\t" + dyr +"\t" + mpx + "\t" + mpy + "\t" + (small+large+xcut+ycut));}
		try {
			target[i+mpx][j+mpy] += large*v;
			if (dx > mpx && dy > mpy)
			{	target[i+mpx+1][j+mpy+1] += small*v; target[i+mpx+1][j+mpy] += xcut*v; target[i+mpx][j+mpy+1] += ycut*v;}
			else if (dx > mpx && dy <= mpy)
			{	target[i+mpx+1][j+mpy-1] += small*v; target[i+mpx+1][j+mpy] += xcut*v; target[i+mpx][j+mpy-1] += ycut*v;}
			else if (dx <= mpx && dy > mpy)
			{	target[i+mpx-1][j+mpy+1] += small*v; target[i+mpx-1][j+mpy] += xcut*v; target[i+mpx][j+mpy+1] += ycut*v;}
			else if (dx <= mpx && dy <= mpy)
			{	target[i+mpx-1][j+mpy-1] += small*v; target[i+mpx-1][j+mpy] += xcut*v; target[i+mpx][j+mpy-1] += ycut*v;}
		}
		catch (ArrayIndexOutOfBoundsException e)
		{
			;
		}
	}
	public static void translateOnePixelAlt(double[][] source, double[][] target, int i, int j, double dx, double dy)
	{
		//first assume that the target position is completely within the range.
		int positiveoffset = 10000;	//this is to make the (int) function round the same way each time.
		double v = source[i][j];
		int N = source.length;
		int M = source[0].length;
		double x = i+dx+0.5, y = j+dy+0.5; //the center of the translated pixel, IF x < 0;
		double cx = (int)(x + positiveoffset) + 0.5, cy = (int)(y + positiveoffset) + 0.5;
		int xmain = (int)cx, ymain = (int)cy;
		cx -= positiveoffset; cy -= positiveoffset;
		xmain -= positiveoffset; ymain -= positiveoffset;
		double deltax = Math.abs(x - cx), deltay = Math.abs(y - cy);
		double big = (1-deltax)*(1-deltay), small = deltax*deltay, xcut = deltax*(1-deltay), ycut = deltay*(1-deltax);
		
		int xoth = x > cx ? xmain+1 : xmain-1;
		int yoth = y > cy ? ymain+1 : ymain-1;
		
//		System.out.print(xmain + "\t" + xoth + "\t" + ymain + "\t" + yoth + "\t" + big + "\t" + xcut);
		//try each one separately.
		if (xmain >= 0 && xmain < N && ymain >= 0 && ymain < M)
			target[xmain][ymain] += v*big;
		if (xoth >= 0 && xoth < N && ymain >= 0 && ymain < M)
			target[xoth][ymain] += v*xcut;
		if (xmain >= 0 && xmain < N && yoth >= 0 && yoth < M)
			target[xmain][yoth] += v*ycut;
		if (xoth >= 0 && xoth < N && yoth >= 0 && yoth < M)
			target[xoth][yoth] += v*small;
	}
	public static void translateOnePixel(double[][] source, double[][] target, int i, int j, int xOff, int yOff, double dx, double dy)
	{
		//This translates pixels to the position plus the offsets.
		//first assume that the target position is completely within the range.
		int N = source.length;
		double v = source[i][j];
		int mpx = dx > 0 ? (int)(dx + 0.5) : (int)(dx - 0.5);
		int mpy = dy > 0 ? (int)(dy + 0.5) : (int)(dy - 0.5);
//		if (i+mpx+xOff <= 1 || i+mpx+xOff >= N-2 || j+mpy+yOff <= 1 || j+mpy+yOff >= N-2) return; //do not translate to the edge.
		double dxr = Math.abs(dx - mpx), dyr = Math.abs(dy - mpy);
		double sdx = 1-dxr, sdy = 1-dyr; 
		double small = dxr*dyr;
		double large = sdx*sdy;
		double xcut = dxr*sdy;
		double ycut = dyr*sdx;
		try {
			target[i+mpx+xOff][j+mpy+yOff] += large*v;
		if (dx > mpx && dy > mpy)
		{	target[i+mpx+1+xOff][j+mpy+1+yOff] += small*v; target[i+mpx+1+xOff][j+mpy+yOff] += xcut*v; target[i+mpx+xOff][j+mpy+1+yOff] += ycut*v;}
		else if (dx > mpx && dy <= mpy)
		{	target[i+mpx+1+xOff][j+mpy-1+yOff] += small*v; target[i+mpx+1+xOff][j+mpy+yOff] += xcut*v; target[i+mpx+xOff][j+mpy-1+yOff] += ycut*v;}
		else if (dx <= mpx && dy > mpy)
		{	target[i+mpx-1+xOff][j+mpy+1+yOff] += small*v; target[i+mpx-1+xOff][j+mpy+yOff] += xcut*v; target[i+mpx+xOff][j+mpy+1+yOff] += ycut*v;}
		else if (dx <= mpx && dy <= mpy)
		{	target[i+mpx-1+xOff][j+mpy-1+yOff] += small*v; target[i+mpx-1+xOff][j+mpy+yOff] += xcut*v; target[i+mpx+xOff][j+mpy-1+yOff] += ycut*v;}
		else
			System.out.println("else");
		}
		catch (ArrayIndexOutOfBoundsException e)
		{
			;
		}
	}
	
	
	//translates a map by the field u which is NxNx2;
	public static double[][] translateMap(double[][] source, double[][][] u)
	{
		int N = source.length; //square;
		int margin = 4, M = N + 2*margin;
		double[][] answer = new double [N][N];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				translateOnePixelAlt(source, answer, i, j, u[i][j][0], u[i][j][1]);
		return answer;
	}
	/**
	 * FFTZ is either null (calculated) or given as fHat2() of the source FFT. Note that t his method modifies FFTZ.
	 * 
	 * @param source
	 * @param fftz
	 * @param r
	 * @return
	 */
	
	public static double[][] translateFancy(double[][] source, double[][][] fftz, double[] r, boolean periodic)
	{
		if (fftz == null)
			fftz = FFTOps.obtainFFT(source).fHat2();
		int nx = fftz.length, ny = fftz[0].length;
		double[][] phase = new double [nx][ny];
		for (int i = 0; i < nx; i++){
			double qx = (i - nx/2)*2*Math.PI/nx;
			for (int j = 0; j < ny; j++){
				double qy = (j - ny/2)*2*Math.PI/ny;
				phase[i][j] = qx*r[0] + qy*r[1];
				fftz[i][j] = Complex.product(fftz[i][j], Complex.expmi(phase[i][j]));
			}
		}
//		LayerViewer.show(Layer.getFreeLayer(phase), 1024, true);
		double[][] ans = new double [nx][ny];
		FFTOps.obtainIFFTCent(fftz, ans);
		if (!periodic){
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
					if (i - r[0] < 0 || i - r[0] > nx-1 || j - r[1] < 0 || j - r[1] > ny-1)
						ans[i][j] = 0;
			FieldOps.changeZeroToAverage(ans);
		}
		return ans;
	}
	/**
	 * This returns the position vector which should make the FFT peaks have the correct phase.
	 * @param source
	 * @param fftz
	 * @return
	 */
	public static double[] getRForTopo(double[][] source, double[][][] fftz, int[][] bragg)
	{
		if (fftz == null)
			fftz = FFTOps.obtainFFT(source).fHat2();
		int nx = fftz.length, ny = fftz[0].length;
		double[] q0 = new double[] {2*Math.PI*bragg[0][0]/nx, 2*Math.PI*bragg[0][1]/ny};
		double[] q1 = new double[] {2*Math.PI*bragg[1][0]/nx, 2*Math.PI*bragg[1][1]/ny};
		double[][] q = new double [][] {q0, q1};
		double phaseq1 = Complex.phaseCenteredZero(fftz[(bragg[0][0] + nx/2) % nx][(bragg[0][1] + ny/2) % ny]);
		double phaseq2 = Complex.phaseCenteredZero(fftz[(bragg[1][0] + nx/2) % nx][(bragg[1][1] + ny/2) % ny]);
	
		double[] phase = new double[] {phaseq1, phaseq2};
		System.out.println("Phase:" + Printer.arrayLnHorizontal(phase));
		RealMatrix ac = new Array2DRowRealMatrix(q, false);
		DecompositionSolver solver = new LUDecomposition(ac).getSolver();
		RealVector bs = new ArrayRealVector(phase, false);
		RealVector solution = null;
		try{solution = solver.solve(bs);
		}
		catch (Exception e)
		{
			System.out.println("Sorry!");
			return null;
		}
		double[] ans = solution.toArray();
		ans[0] = ans[0];
		ans[1] = ans[1];
		System.out.println("" + ans[0] + " * " + q0[0] + " = " + ans[0]*q0[0]);
		System.out.println("" + ans[1] + " * " + q0[1] + " = " + ans[1]*q0[1]);
		System.out.println("" + ans[0] + " * " + q0[0] + " + " + ans[1] + " * " + q0[1] + " = " + (ans[0]*q0[0] + ans[1]*q0[1]));
		System.out.println("" + ans[0] + " * " + q1[0] + " = " + ans[0]*q1[0]);
		System.out.println("" + ans[1] + " * " + q1[1] + " = " + ans[1]*q1[1]);
		System.out.println("" + ans[0] + " * " + q1[0] + " + " + ans[1] + " * " + q1[1] + " = " + (ans[0]*q1[0] + ans[1]*q1[1]));
		
		return ans;
	}
//		int N = source.length; //square;
//		int margin = 4, M = N + 2*margin;
//		double[][] answer = new double [N][N];
//		double[][] A = {
//				{sxx, sxy, sx},
//				{sxy, syy, sy},
//				{sx, sy, n}};
//		double[] B = {sxz,syz,sz};		
		//use Apahce linear solver:
//		RealMatrix ac = new Array2DRowRealMatrix(A, false);
//		DecompositionSolver solver = new LUDecomposition(ac).getSolver();
//		RealVector bs = new ArrayRealVector(B, false);
	public static double[][] applyUField(double[][] source, double[][][] u, int detail, int nsec)
	{
		int N = source.length, M = source[0].length;
		double[][] target = new double [N][M];
		int secsize = N/nsec;
		double[][] temp = new double [N*detail/nsec][M*detail/nsec];
		double mean = FieldOps.mean(source);
		double[][] tempsec = new double [N/nsec][M/nsec];
		double a, b; int ai, bi;
		double[] r = new double [2];
		double dr = 1/(double)(detail), edge = 1/(double)(2*detail);
		for (int m = 0; m < nsec; m++){
			for (int n = 0; n < nsec; n++){
				for (int i = 0; i < N*detail/nsec; i++)
					for (int j = 0; j < M*detail/nsec; j++)
					{
						r[0] = edge + i*dr + m*secsize;
						r[1] = edge + j*dr + n*secsize;
						a = r[0] - u[(int)r[0]][(int)r[1]][0];
						b = r[1] - u[(int)r[0]][(int)r[1]][1];
//						ai = roundDown(a); bi = roundDown(b);
//						if (ai >= 0 && ai < N && bi >= 0 && bi < M)
//							temp[i][j] = source[ai][bi];
//						else
//							temp[i][j] = 0;
						temp[i][j] = getValueAt(source, a, b, mean);
					}
				reduce(detail, temp, tempsec);
				for (int i = 0; i < tempsec.length; i++)
					for (int j = 0; j < tempsec[0].length; j++)
						target[m*secsize + i][n*secsize + j] = tempsec[i][j];
			}
//			System.out.print(m + " ");
		}
		return target;
	}
	public static double[][] applyUFieldSmooth(double[][] source, double[][][] u, int detail, int nsec)
	{
		int N = source.length, M = source[0].length;
		double[][] target = new double [N][M];
		int secsize = N/nsec;
		double[][] temp = new double [N*detail/nsec][M*detail/nsec];
		double[][] tempsec = new double [N/nsec][M/nsec];
		double a, b; int ai, bi;
		double mean = FieldOps.mean(source);
		double[] r = new double [2];
		double dr = 1/(double)(detail), edge = 1/(double)(2*detail);
		for (int m = 0; m < nsec; m++){
			for (int n = 0; n < nsec; n++){
				for (int i = 0; i < N*detail/nsec; i++)
					for (int j = 0; j < M*detail/nsec; j++)
					{
						r[0] = edge + i*dr + m*secsize;
						r[1] = edge + j*dr + n*secsize;
						a = r[0] - getValueAt(u, r[0], r[1], 0);
						b = r[1] - getValueAt(u, r[0], r[1], 1);
//						ai = roundDown(a); bi = roundDown(b);
//						if (ai >= 0 && ai < N && bi >= 0 && bi < M)
//							temp[i][j] = getValueAt(source, a, b);
//						else
//							temp[i][j] = 0;
						temp[i][j] = getValueAt(source, a, b, mean);
					}
				reduce(detail, temp, tempsec);
				for (int i = 0; i < tempsec.length; i++)
					for (int j = 0; j < tempsec[0].length; j++)
						target[m*secsize + i][n*secsize + j] = tempsec[i][j];
			}
//			System.out.print(m + " ");
		}
		return target;
	}
	public static void applyUFieldSmooth(double[][] source, double[][][] u, int detail, int nsec, double[][] target)
	{
		int N = source.length, M = source[0].length;
		int secsize = N/nsec;
		double[][] temp = new double [N*detail/nsec][M*detail/nsec];
		double[][] tempsec = new double [N/nsec][M/nsec];
		double a, b; int ai, bi;
		double[] r = new double [2];
//		double mean = 0;//FieldOps.mean(source);
		double mean = FieldOps.mean(source);
		
		double dr = 1/(double)(detail), edge = 1/(double)(2*detail);
		for (int m = 0; m < nsec; m++){
			for (int n = 0; n < nsec; n++){
				for (int i = 0; i < N*detail/nsec; i++)
					for (int j = 0; j < M*detail/nsec; j++)
					{
						r[0] = edge + i*dr + m*secsize;
						r[1] = edge + j*dr + n*secsize;
						a = r[0] - getValueAt(u, r[0], r[1], 0);
						b = r[1] - getValueAt(u, r[0], r[1], 1);
						ai = roundDown(a); bi = roundDown(b);
//						if (ai >= 0 && ai < N && bi >= 0 && bi < M)
//							temp[i][j] = getValueAt(source, a, b);
//						else
//							temp[i][j] = 0;
						temp[i][j] = getValueAt(source, a, b, mean);
					}
				reduce(detail, temp, tempsec);
				for (int i = 0; i < tempsec.length; i++)
					for (int j = 0; j < tempsec[0].length; j++)
						target[m*secsize + i][n*secsize + j] = tempsec[i][j];
			}
//			System.out.print(m + " ");
		}
	}
	public static void applyUFieldBiCubic(double[][] source, double[][][] u, double[][] target)
	{
		int N = source.length, M = source[0].length;
		double[] r = new double [2];
//		double mean = 0;//FieldOps.mean(source);
		double mean = FieldOps.mean(source);
		double a, b;
		
		
		BicubicSplineInterpolatingFunction interp = FieldOps.getBicubicInterpoation(source);
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
			{
				r[0] = i;
				r[1] = j;
				a = r[0] - u[i][j][0];
				b = r[1] - u[i][j][1];
				if (a >= 0 && a <= N-1 && b >= 0 && b <= M-1)
					target[i][j] = interp.value(a, b);
				else
					target[i][j] = mean;
			}
//			System.out.print(m + " ");
		}
	/**
	 * Uses a pre-initialized interpolator
	 * @param interp
	 * @param u
	 * @param target
	 */
	public static void applyUFieldBiCubic(BicubicSplineInterpolatingFunction interp, double[][][] u, double[][] target, double mean)
	{
		int N = target.length, M = target[0].length;
		double[] r = new double [2];
//		double mean = 0;//FieldOps.mean(source);
		double a, b;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
			{
				r[0] = i;
				r[1] = j;
				a = r[0] - u[i][j][0];
				b = r[1] - u[i][j][1];
				if (a >= 0 && a <= N-1 && b >= 0 && b <= M-1)
					target[i][j] = interp.value(a, b);
				else
					target[i][j] = mean;
			}
//			System.out.print(m + " ");
		}
	public static boolean[][] getOutsidePixels(double[][][] u)
	{
		int N = u.length, M = u[0].length;
		double[] r = new double [2];
//		double mean = 0;//FieldOps.mean(source);
		double a, b;
		boolean[][] outside = new boolean [N][M];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
			{
				r[0] = i;
				r[1] = j;
				a = r[0] - u[i][j][0];
				b = r[1] - u[i][j][1];
					outside[i][j] = !(a >= 0 && a < N && b >= 0 && b < M);
			}
		return outside;
//			System.out.print(m + " ");
		}

	//translates all pixles between i and i+di, j and j+dj by (dx, dy);
	public static void translatePixelGroup(double[][] source, double[][] target, int i, int j, int di, int dj, double dx, double dy)
	{
//		System.out.print("(" + dx + ", " + dy + ")");
		for (int m = i; m < i+di; m++)
			for (int n = j; n < j+dj; n++)
			{
				translateOnePixel(source, target, m, n, dx, dy);
			}
	}
	public static double[][] reduce(int factor, double[][] data)
	{
		double[][] field = new double [data.length/factor][data[0].length/factor];
		for (int i = 0; i < field.length; i++)
		for(int j = 0; j < field[0].length; j++)
		{
			field[i][j] = 0;
			for (int k = 0; k < factor; k++)
				for (int m = 0; m < factor; m++)
				{
					field[i][j] += data[i*factor + k][j*factor + m];
				}
			field[i][j] /= factor*factor;
		}
		return field;

	}
	
	/**
	 * this assumes data of the form [k][i][j] NOT [i][j][k]. Note order of the parameters lk, li, lj.
	 * The convention is due to the fact that in 3D maps the data is stored [nlayers][nx][ny] and i and j are x and y.
	 * @param data
	 * @param li
	 * @param lj
	 * @param lk
	 * @return
	 */
	public static double[][][] getGaussianSmoothing3D(double[][][] data, double lk, double li, double lj)
	{
		int gli = (int)(3*(li+1));
		int glj = (int)(3*(lj+1));
		int glk = (int)(3*(lk+1));
		int nx = data[0].length;
		int ny = data[0][0].length;
		int nlayers = data.length;
		double[][][] smooth = new double[nlayers][nx][ny];
		
		if (li <= 0 || lj <= 0 && lk > 0){
			double[] tempspec = new double [nlayers];
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
				{
					for (int k = 0; k < nlayers; k++)
						tempspec[k] = data[k][i][j];
					tempspec = ArrayOps.gaussSmooth(tempspec, lk);
					for (int k = 0; k < nlayers; k++)
						smooth[k][i][j] = tempspec[k];
				}
			return smooth;
		}
		else if (lk <= 0)
			return data.clone();
		
		//now assuming li and lj are both positive
		if (lk <= 0){
			for (int k = 0; k < nlayers; k++)
				smooth[k] = FieldOps.gaussSmooth(data[k], li, lj);
			return smooth;
		}
		
		double[][][] gauss = new double [glk+1][gli+1][glj+1];

		for (int k = 0; k < glk+1; k++)
			for (int i = 0; i < gli+1; i++)
				for (int j = 0; j < glj+1; j++)
					gauss[k][i][j] = Math.exp(-((k*k)/(2*lk*lk) + (i*i)/(2*li*li) + (j*j)/(2*lj*lj)));
		
		double gsum;
		int xprime, yprime, zprime, xpmin, xpmax, ypmin, ypmax, zpmin, zpmax;
		for (int k = 0; k < nlayers; k++){System.out.print(" " + k);
			for (int i = 0; i < nx; i++){
				for (int j = 0; j < ny; j++)
				{
					gsum = 0;
					smooth[k][i][j] = 0;
					//this is done so that xprime runs from 0 to N, or from x - biggerL/2 to x+biggerL/2
					xpmin = Math.max(0, i - gli);
					xpmax = Math.min(nx, i + gli);
					ypmin = Math.max(0, j - glj);
					ypmax = Math.min(ny, j + glj);
					zpmin = Math.max(0, k - glk);
					zpmax = Math.min(nlayers, k + glk);
					for (zprime = zpmin; zprime < zpmax; zprime++)
						for (xprime = xpmin; xprime < xpmax; xprime++)
							for (yprime = ypmin; yprime < ypmax; yprime++)
							{
								gsum += gauss[Math.abs(k-zprime)][Math.abs(i-xprime)][Math.abs(j-yprime)];
								smooth[k][i][j] += gauss[Math.abs(k-zprime)][Math.abs(i-xprime)][Math.abs(j-yprime)]*data[zprime][xprime][yprime];
							}
					smooth[k][i][j] /= gsum;
				}
			}
		}
		return smooth;
	}
	public static void rescale(double[][] data, double min, double max)
	{
		double dmin = ArrayOps.min(data), dmax = ArrayOps.max(data);
		double norm;
		double ddata = dmax - dmin, delta = max - min;
		for (int i = 0; i < data.length; i++)
			for(int j = 0; j < data[0].length; j++)
			{
				norm = (data[i][j] - dmin)/ddata;
				data[i][j] = min + norm*delta;
			}
	}
	public static double[][] rescaleXY(double[][] data, double scaleFactor, int detail)
	{
		double[][] ans = new double[data.length][data[0].length];
		FieldOps.applyLinTransNum(data, ans, new double[][] {{scaleFactor, 0}, {0, scaleFactor}}, new int[] {data.length/2,data[0].length/2}, detail);
		return ans;
	}
	public static double paraSqOver(int x, int y, double a, double b, double[] aHat, double[] bHat)
	{
		double a00 = DataManip.dot(aHat, x, y);
		double a01 = DataManip.dot(aHat, x+1, y);
		double a10 = DataManip.dot(aHat, x, y+1);
		double a11 = DataManip.dot(aHat, x+1, y+1);
		double b00 = DataManip.dot(bHat, x, y);
		double b01 = DataManip.dot(bHat, x+1, y);
		double b10 = DataManip.dot(bHat, x, y+1);
		double b11 = DataManip.dot(bHat, x+1, y+1);
		
		int ai00 = roundDown(a00);
		int ai10 = roundDown(a10);
		int ai01 = roundDown(a01);
		int ai11 = roundDown(a11);
		int bi00 = roundDown(b00);
		int bi10 = roundDown(b10);
		int bi01 = roundDown(b01);
		int bi11 = roundDown(b11);
		
		return 0;
	}
	
	//applies the linear transformation numerically, as exactly as specified
	public static void applyLinTransNum(double[][] source, double[][] target, double[][] matrix, int[] origin, int detail)
	{
		int N = source.length, M = source[0].length;
		double[][] temp = new double [N*detail][M*detail];
		double[][] inv = MatrixOps.inv2x2(matrix);
		double[] aHat = inv[0];//new double[] {inv[0][0], inv[1][0]};
		double[] bHat = inv[1];//new double[] {inv[0][1], inv[1][1]};
		double a, b; int ai, bi;
		double[] r = new double [2];
		double dr = 1/(double)(detail), edge = 1/(double)(2*detail);
		for (int i = 0; i < N*detail; i++)
			for (int j = 0; j < M*detail; j++)
			{
				r[0] = edge + i*dr - origin[0];
				r[1] = edge + j*dr - origin[1];
				a = DataManip.dot(r, aHat);
				b = DataManip.dot(r, bHat);
				ai = roundDown(a)+origin[0]; bi = roundDown(b)+origin[1];
				if (ai >= 0 && ai < N && bi >= 0 && bi < M)
					temp[i][j] = source[ai][bi];
				else
					temp[i][j] = 0;
			}
		reduce(detail, temp, target);
	}

	//applies the linear transformation in sections, so no heap overflow:
	public static void applyLinTransNumSec(double[][] source, double[][] target, double[][] matrix, int[] origin, int detail, int nsec)
	{
		int N = source.length, M = source[0].length;
		int Nt = target.length, Mt = target[0].length;
		int secsize = Nt/nsec;
		double[][] temp = new double [Nt*detail/nsec][Mt*detail/nsec];
		double[][] tempsec = new double [Nt/nsec][Mt/nsec];
		double[][] inv = MatrixOps.inv2x2(matrix);
		double[] aHat = inv[0];//new double[] {inv[0][0], inv[1][0]};
		double[] bHat = inv[1];//new double[] {inv[0][1], inv[1][1]};
		double a, b; int ai, bi;
		double[] r = new double [2];
		double dr = 1/(double)(detail), edge = 1/(double)(2*detail);
		for (int m = 0; m < nsec; m++){
			for (int n = 0; n < nsec; n++){
				for (int i = 0; i < Nt*detail/nsec; i++)
					for (int j = 0; j < Mt*detail/nsec; j++)
					{
						r[0] = edge + i*dr + m*secsize - origin[0];
						r[1] = edge + j*dr + n*secsize - origin[1];
						a = DataManip.dot(r, aHat);
						b = DataManip.dot(r, bHat);
						ai = roundDown(a)+origin[0]; bi = roundDown(b)+origin[1];
						if (ai >= 0 && ai < N && bi >= 0 && bi < M)
							temp[i][j] = source[ai][bi];
						else
							temp[i][j] = 0;
					}
				reduce(detail, temp, tempsec);
				for (int i = 0; i < tempsec.length; i++)
					for (int j = 0; j < tempsec[0].length; j++)
						target[m*secsize + i][n*secsize + j] = tempsec[i][j];
			}
//			System.out.print(m + " ");
		}
	}
	public static void applyLinTransNumSecBi(double[][] source, double[][] target, double[][] matrix, int[] origin, int detail, int nsec)
	{
		int N = source.length, M = source[0].length;
		int Nt = target.length, Mt = target[0].length;
		int secsize = Nt/nsec;
		double[][] temp = new double [Nt*detail/nsec][Mt*detail/nsec];
		double[][] tempsec = new double [Nt/nsec][Mt/nsec];
		double[][] inv = MatrixOps.inv2x2(matrix);
		double[] aHat = inv[0];//new double[] {inv[0][0], inv[1][0]};
		double[] bHat = inv[1];//new double[] {inv[0][1], inv[1][1]};
		double a, b; int ai, bi;
		double mean = FieldOps.mean(source);
		double[] r = new double [2];
		double dr = 1/(double)(detail), edge = 1/(double)(2*detail);
		for (int m = 0; m < nsec; m++){
			for (int n = 0; n < nsec; n++){
				for (int i = 0; i < Nt*detail/nsec; i++)
					for (int j = 0; j < Mt*detail/nsec; j++)
					{
						r[0] = edge + i*dr + m*secsize - origin[0];
						r[1] = edge + j*dr + n*secsize - origin[1];
						a = DataManip.dot(r, aHat);
						b = DataManip.dot(r, bHat);
						ai = roundDown(a)+origin[0]; bi = roundDown(b)+origin[1];
						if (ai >= 0 && ai < N && bi >= 0 && bi < M)
							temp[i][j] = getValueAt(source, a, b, mean);
						else
							temp[i][j] = 0;
					}
				reduce(detail, temp, tempsec);
				for (int i = 0; i < tempsec.length; i++)
					for (int j = 0; j < tempsec[0].length; j++)
						target[m*secsize + i][n*secsize + j] = tempsec[i][j];
			}
//			System.out.print(m + " ");
		}
	}
	public static void applyLinTransNumSec(double[][] source, double[][] target, double[][] matrix, int[] origin, int detail, int nsec, int[] targorigin)
	{
		int N = source.length, M = source[0].length;
		int Nt = target.length, Mt = target[0].length;
		int secsize = Nt/nsec;
		double[][] temp = new double [Nt*detail/nsec][Mt*detail/nsec];
		double[][] tempsec = new double [Nt/nsec][Mt/nsec];
		double[][] inv = MatrixOps.inv2x2(matrix);
		double[] aHat = inv[0];//new double[] {inv[0][0], inv[1][0]};
		double[] bHat = inv[1];//new double[] {inv[0][1], inv[1][1]};
		double a, b; int ai, bi;
		double[] r = new double [2];
		double dr = 1/(double)(detail), edge = 1/(double)(2*detail);
		for (int m = 0; m < nsec; m++){
			for (int n = 0; n < nsec; n++){
				for (int i = 0; i < Nt*detail/nsec; i++)
					for (int j = 0; j < Mt*detail/nsec; j++)
					{
						r[0] = edge + i*dr + m*secsize - targorigin[0];
						r[1] = edge + j*dr + n*secsize - targorigin[1];
						a = DataManip.dot(r, aHat);
						b = DataManip.dot(r, bHat);
						ai = roundDown(a)+origin[0]; bi = roundDown(b)+origin[1];
						if (ai >= 0 && ai < N && bi >= 0 && bi < M)
							temp[i][j] = source[ai][bi];
						else
							temp[i][j] = 0;
					}
				reduce(detail, temp, tempsec);
				for (int i = 0; i < tempsec.length; i++)
					for (int j = 0; j < tempsec[0].length; j++)
						target[m*secsize + i][n*secsize + j] = tempsec[i][j];
			}
//			System.out.print(m + " ");
		}
	}
	public static void applyLinTransNumBiggerSec(double[][] source, double[][] target, double[][] matrix, int[] origin, int detail, int nsec)
	{
		int bigfactor = target.length/source.length;
		//the source needs to be blown up to meet the target:
		double[][] bigmat = new double [2][2];
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				bigmat[i][j] = matrix[i][j]*bigfactor;
		applyLinTransNumSec(source, target, bigmat, origin, detail, nsec, new int [] {origin[0]*bigfactor, origin[1]*bigfactor});
	}
	
	//This one first expands the source array, then performs the transformation at some detail, then contracts again.
	public static void applyLinTransNum2Sec(double[][] source, double[][] target, double[][] matrix, int[] origin, int detail, int nsec, int expscale, int nexpsec)
	{
		int N = source.length;
		int Nxp = N*expscale/nexpsec;
		//the expansion is in a power of 2, the number of sections to be expanded;
		double[][] tempexp = null;
		double[][] temptarg = null;
		double[][] temp = null;
		boolean exp2 = expscale == 2;
		switch(expscale)
		{
		case 1:
			FieldOps.applyLinTransNumSec(source, target, matrix, origin, detail, nsec);
			return;
		case 2: case 4:
//			tempexp = new double [Nxp][Nxp];
			temptarg = new double [Nxp][Nxp];
			temp = new double [N/nexpsec][N/nexpsec];
			break;
		default:
			return;
		}
		
		int dx = N/nexpsec, dy = N/nexpsec;
		//this is the origin in the coordinate system of tempexp.
		int[] neworigin = new int [2];
//		int[] exporigin = {origin[0]*expscale, origin[1]*expscale};
		for (int i = 0; i < nexpsec; i++)
			for (int j = 0; j < nexpsec; j++)
			{
				tempexp = exp2 ? FieldOps.expand(source, i*dx, (i+1)*dx, j*dy, (j+1)*dy) : FieldOps.expand4(source, i*dx, (i+1)*dx, j*dy, (j+1)*dy);
				neworigin[0] = (origin[0]-i*dx)*expscale;
				neworigin[1] = (origin[1]-j*dy)*expscale;
				FieldOps.applyLinTransNumSec(tempexp, temptarg, matrix, neworigin, detail, nsec);
				FieldOps.reduce(expscale, temptarg, temp);
				for (int m = 0; m < dx; m++)
					for (int n = 0; n < dy; n++)
						target[i*dx+m][j*dy+n] = temp[m][n];
			}
	}
	public static void applyLinTransNumSecComplex(double[][][] source, double[][][] target, double[][] matrix, int[] origin, int detail, int nsec)
	{
		int N = source.length, M = source[0].length;
		int secsize = N/nsec;
		double[][][] temp = new double [N*detail/nsec][M*detail/nsec][2];
		double[][][] tempsec = new double [N/nsec][M/nsec][2];
		double[][] inv = MatrixOps.inv2x2(matrix);
		double[] aHat = inv[0];//new double[] {inv[0][0], inv[1][0]};
		double[] bHat = inv[1];//new double[] {inv[0][1], inv[1][1]};
		double a, b; int ai, bi;
		double[] r = new double [2];
		double dr = 1/(double)(detail), edge = 1/(double)(2*detail);
		for (int m = 0; m < nsec; m++){
			for (int n = 0; n < nsec; n++){
				for (int i = 0; i < N*detail/nsec; i++)
					for (int j = 0; j < M*detail/nsec; j++)
					{
						r[0] = edge + i*dr + m*secsize - origin[0];
						r[1] = edge + j*dr + n*secsize - origin[1];
						a = DataManip.dot(r, aHat);
						b = DataManip.dot(r, bHat);
						ai = roundDown(a)+origin[0]; bi = roundDown(b)+origin[1];
						if (ai >= 0 && ai < N && bi >= 0 && bi < M)
							{temp[i][j][0] = source[ai][bi][0]; temp[i][j][1] = source[ai][bi][1];}
						else
							{temp[i][j][0] = 0; temp[i][j][1] = 0;}
					}
				reduce(detail, temp, tempsec);
				for (int i = 0; i < tempsec.length; i++)
					for (int j = 0; j < tempsec[0].length; j++)
						{target[m*secsize + i][n*secsize + j][0] = tempsec[i][j][0]; target[m*secsize + i][n*secsize + j][1] = tempsec[i][j][1];}
			}
//			System.out.print(m + " ");
		}
	}

	
	public static void applyLinTransNum(double[][] source, double[][] target, double[][] matrix, int[] origin, double[][] temp)
	{
		int N = source.length, M = source[0].length;
		int detail = temp.length/N;
		double[][] inv = MatrixOps.inv2x2(matrix);
		double[] aHat = inv[0];//new double[] {inv[0][0], inv[1][0]};
		double[] bHat = inv[1];//new double[] {inv[0][1], inv[1][1]};
		double a, b; int ai, bi;
		double[] r = new double [2];
		double dr = 1/(double)(detail), edge = 1/(double)(2*detail);
		for (int i = 0; i < N*detail; i++)
			for (int j = 0; j < M*detail; j++)
			{
				r[0] = edge + i*dr - origin[0];
				r[1] = edge + j*dr - origin[1];
				a = DataManip.dot(r, aHat);
				b = DataManip.dot(r, bHat);
				ai = roundDown(a)+origin[0]; bi = roundDown(b)+origin[1];
				if (ai >= 0 && ai < N && bi >= 0 && bi < M)
					temp[i][j] = source[ai][bi];
				else
					temp[i][j] = 0;
			}
		reduce(detail, temp, target);
	}
	
	public static int roundDown(double x)
	{
		return (int)(x+1000000)-1000000;
	}
	public static int round(double x)
	{
		return (int)(x+1000000+0.5)-1000000;
	}
	//rounds to the nearest even number
	public static int roundEven(double x)
	{
		return round(x/2)*2;
	}
	//returns atan(y/x) as an angle between 0 and 2pi
	public static void reduce(int factor, double[][] data, double[][] target)
	{
		for (int i = 0; i < target.length; i++)
		for(int j = 0; j < target[0].length; j++)
		{
			target[i][j] = 0;
			for (int k = 0; k < factor; k++)
				for (int m = 0; m < factor; m++)
				{
					target[i][j] += data[i*factor + k][j*factor + m];
				}
			target[i][j] /= factor*factor;
		}
	}
	public static void reduce(int factor, double[][][] data, double[][][] target)
	{
		for (int i = 0; i < target.length; i++)
		for(int j = 0; j < target[0].length; j++)
		{
			target[i][j][0] = 0;
			target[i][j][1] = 0;
			for (int k = 0; k < factor; k++)
				for (int m = 0; m < factor; m++)
				{
					target[i][j][0] += data[i*factor + k][j*factor + m][0];
					target[i][j][1] += data[i*factor + k][j*factor + m][1];
					
				}
			target[i][j][0] /= factor*factor;
			target[i][j][1] /= factor*factor;
		}
	}
	public static double[][][] reduce(int factor, double[][][] data)
	{
		double[][][] target = new double [data.length/factor][data[0].length/factor][2];
		for (int i = 0; i < target.length; i++)
		for(int j = 0; j < target[0].length; j++)
		{
			target[i][j][0] = 0;
			target[i][j][1] = 0;
			for (int k = 0; k < factor; k++)
				for (int m = 0; m < factor; m++)
				{
					target[i][j][0] += data[i*factor + k][j*factor + m][0];
					target[i][j][1] += data[i*factor + k][j*factor + m][1];
					
				}
			target[i][j][0] /= factor*factor;
			target[i][j][1] /= factor*factor;
		}
		return target;
	}

	public static void translatePixelGroup(double[][] source, double[][] target, int i, int j, int di, int dj, int xOff, int yOff, double dx, double dy)
	{
//		System.out.print("(" + dx + ", " + dy + ")");
		for (int m = i; m < i+di; m++)
			for (int n = j; n < j+dj; n++)
			{
				translateOnePixel(source, target, m, n, xOff, yOff, dx, dy);
			}
	}
	
	
	//This method translates square blocks.
	public static double[][] translateMapBlocks(double[][] topo, double[][][] u, int blocksize)
	{
		int N = topo.length;
		int center = blocksize/2;
		int n = N/blocksize;
		int margin = 4;
		double[][] target = new double [N+2*margin][N+2*margin];
		double[][] topocorr = new double [N][N];
		int x, y;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
			{
				x = i*blocksize;
				y = j*blocksize;
				translatePixelGroup(topo, target, x, y, blocksize, blocksize, margin, margin, u[x+center][y+center][0], u[x+center][y+center][1]);
			}
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				topocorr[i][j] = target[i+margin][j+margin];
		return topocorr;
	}
	
	//frequency must be a power of 2.
	public static double[][] gaussSmooth(double[][] data, int L, int frequency)
	{
		int biggerL = 4*L;
		int N = data.length;
		double[][] gauss = new double [biggerL][biggerL];
		int size = N/frequency;
		double[][] smooth = new double[size][size];
		for (int k = 0; k < biggerL; k++)
			for (int l = 0; l < biggerL; l++)
				gauss[k][l] = Math.exp((-((double)(k*k) + (l*l))/(2*L*L)));
		
		double gsum;
		int x, y;
		int xprime, yprime, xpmin, xpmax, ypmin, ypmax;
		for (int i = 0; i < size; i++){if (i % 10 == 0) System.out.println(i + " out of " + size);
			for (int j = 0; j < size; j++)
			{
				gsum = 0;
				smooth[i][j] = 0;
				x = i*frequency;
				y = j*frequency;
				//this is done so that xprime runs from 0 to N, or from x - biggerL/2 to x+biggerL/2
				xpmin = Math.max(0, x - 3*L);
				xpmax = Math.min(N, x + 3*L);
				ypmin = Math.max(0, y - 3*L);
				ypmax = Math.min(N, y + 3*L);
				for (xprime = xpmin; xprime < xpmax; xprime++)
					for (yprime = ypmin; yprime < ypmax; yprime++)
					{
						gsum += gauss[Math.abs(x-xprime)][Math.abs(y-yprime)];
						smooth[i][j] += gauss[Math.abs(x-xprime)][Math.abs(y-yprime)]*data[xprime][yprime];
					}
				smooth[i][j] /= gsum;
			}
		}
		while(size < N){
			smooth = FieldOps.expand(smooth);
			size *= 2;
		}
		return smooth;
	}
	public static double[][] gaussSmooth(double[][] data, double L)
	{
		int gl = (int)(3*(L+1));
		int N = data.length;
		int M = data[0].length;
		double[][] smooth = new double[N][M];
		
		if (L == 0){
			return data.clone();
		}
		if (L == Double.POSITIVE_INFINITY || L == -1)
		{
			double d = mean(data);
			for (int i = 0; i < N; i++)//if (i % 10 == 0) System.out.println(i + " out of " + size);
				for (int j = 0; j < M; j++)
					smooth[i][j] = d;
			return smooth;
		}
		
		double[][] gauss = new double [gl+1][gl+1];

		for (int k = 0; k < gl+1; k++)
			for (int l = 0; l < gl+1; l++)
				gauss[k][l] = Math.exp(-((k*k) + (l*l))/(2*L*L));
		
		double gsum;
		int xprime, yprime, xpmin, xpmax, ypmin, ypmax;
		for (int i = 0; i < N; i++){//if (i % 10 == 0) System.out.println(i + " out of " + size);
			for (int j = 0; j < M; j++)
			{
				gsum = 0;
				smooth[i][j] = 0;
				//this is done so that xprime runs from 0 to N, or from x - biggerL/2 to x+biggerL/2
				xpmin = Math.max(0, i - gl);
				xpmax = Math.min(N, i + gl);
				ypmin = Math.max(0, j - gl);
				ypmax = Math.min(M, j + gl);
				for (xprime = xpmin; xprime < xpmax; xprime++)
					for (yprime = ypmin; yprime < ypmax; yprime++)
					{
						gsum += gauss[Math.abs(i-xprime)][Math.abs(j-yprime)];
						smooth[i][j] += gauss[Math.abs(i-xprime)][Math.abs(j-yprime)]*data[xprime][yprime];
					}
				smooth[i][j] /= gsum;
			}
		}
		return smooth;
	}
	public static void gaussSmooth(double[][] data, double[][] smooth, double L)
	{
		int gl = (int)(3*(L+1));
		int N = data.length;
		int M = data[0].length;
	
		if (L == Double.POSITIVE_INFINITY || L == -1)
		{
			double d = mean(data);
			for (int i = 0; i < N; i++)//if (i % 10 == 0) System.out.println(i + " out of " + size);
				for (int j = 0; j < M; j++)
					smooth[i][j] = d;
		}
		
		double[][] gauss = new double [gl+1][gl+1];

		for (int k = 0; k < gl+1; k++)
			for (int l = 0; l < gl+1; l++)
				gauss[k][l] = Math.exp(-((k*k) + (l*l))/(2*L*L));
		
		double gsum;
		int xprime, yprime, xpmin, xpmax, ypmin, ypmax;
		for (int i = 0; i < N; i++){//if (i % 10 == 0) System.out.println(i + " out of " + size);
			for (int j = 0; j < M; j++)
			{
				gsum = 0;
				smooth[i][j] = 0;
				//this is done so that xprime runs from 0 to N, or from x - biggerL/2 to x+biggerL/2
				xpmin = Math.max(0, i - gl);
				xpmax = Math.min(N, i + gl);
				ypmin = Math.max(0, j - gl);
				ypmax = Math.min(M, j + gl);
				for (xprime = xpmin; xprime < xpmax; xprime++)
					for (yprime = ypmin; yprime < ypmax; yprime++)
					{
						gsum += gauss[Math.abs(i-xprime)][Math.abs(j-yprime)];
						smooth[i][j] += gauss[Math.abs(i-xprime)][Math.abs(j-yprime)]*data[xprime][yprime];
					}
				smooth[i][j] /= gsum;
			}
		}
	}
	public static double[][] getGaussSmoothConvolve(double[][] data, double L)
	{
		int N = data.length;
		int M = data[0].length;
		double[][] smooth = new double[N][M];
		
		if (L == 0){
			smooth[0][0] = 1;
			return smooth;
		}
		if (L == Double.POSITIVE_INFINITY || L == -1)
		{
			for (int i = 0; i < N; i++)//if (i % 10 == 0) System.out.println(i + " out of " + size);
				for (int j = 0; j < M; j++)
					smooth[i][j] = 1.0d/(N*M);
			return smooth;
		}
		
		
		for (int k = 0; k < N; k++)
			for (int l = 0; l < M; l++)
				smooth[k][l] = Math.exp(-((k-N/2)*(k-N/2) + (l-M/2)*(l-M/2))/(2*L*L));
		
		double[][] smoothAns = new double[N][M];
		FieldOps.shift(smooth, smoothAns);
		
		double sum = FieldOps.sum(smoothAns);		
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				smoothAns[i][j]/=sum;
		
		return smoothAns;
	}
	/**
	 * This is intended to center the positive lobe on the +/- a axis of the latice, the negative lobe 90 degrees (hopefully the +/- b axis)
	 * The radial function is to be r*exp(r/L)
	 * @param data
	 * @param L
	 * @param latt
	 * @return
	 */
	public static double[][] getDWaveConvolve(double[][] data, double L, AtomicCoordinatesSet latt)
	{
		int gl = (int)(3*(L+1));
		int N = data.length;
		int M = data[0].length;
		
		double[][] smooth = new double [N][M];
		double[][] radial = new double [N][M];
		double lattPhi = latt.getAngleBetween(new double[] {1, 0}, 0);
//		System.out.println(lattPhi);
		for (int k = 0; k < N; k++)
			for (int l = 0; l < M; l++)
			{
				double r = Distance.distance(k-N/2, l-M/2);
				double phi = atan(k-N/2, l-M/2) - lattPhi;
				radial[k][l] = r*Math.exp(-r/L);
				smooth[k][l] = radial[k][l]*Math.cos(2*phi);
			}
		double[][] smoothAns = new double[N][M];
		FieldOps.shift(smooth, smoothAns);
		
		double sum = FieldOps.sum(radial);		
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				smoothAns[i][j]/=sum;
		
		return smoothAns;
	}
	public static double[][] getDWaveConvolve_sharper(double[][] data, double L, AtomicCoordinatesSet latt)
	{
		int gl = (int)(3*(L+1));
		int N = data.length;
		int M = data[0].length;
		
		double[][] smooth = new double [N][M];
		double[][] radial = new double [N][M];
		double lattPhi = latt.getAngleBetween(new double[] {1, 0}, 0);
//		System.out.println(lattPhi);
		for (int k = 0; k < N; k++)
			for (int l = 0; l < M; l++)
			{
				double r = Distance.distance(k-N/2, l-M/2);
				double phi = (atan(k-N/2, l-M/2) - lattPhi + 2*Math.PI) % 2*Math.PI;
				double angular = phi > Math.PI ? (-1 + 2/(1+Math.exp((phi-3*Math.PI/2)/5))) : (1 - 2/(1+Math.exp((phi-Math.PI/2)/5)));
				radial[k][l] = r*Math.exp(-r/L);
				smooth[k][l] = radial[k][l]*angular;
			}
		double[][] smoothAns = new double[N][M];
		FieldOps.shift(smooth, smoothAns);
		
		double sum = FieldOps.sum(radial);		
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				smoothAns[i][j]/=sum;
		
		return smoothAns;
	}
	/**
	 * This is intended to center the nodal line on the +/- a axis of the latice, so as to distinguish one type of swastika from its mirror image.
	 * @param data
	 * @param L
	 * @param latt
	 * @return
	 */
	public static double[][] getGWaveConvolve(double[][] data, double L, AtomicCoordinatesSet latt)
	{
		int gl = (int)(3*(L+1));
		int N = data.length;
		int M = data[0].length;
		
		double[][] smooth = new double [N][M];
		double[][] radial = new double [N][M];
		double lattPhi = latt.getAngleBetween(new double[] {1, 0}, 0);
//		System.out.println(lattPhi);
		for (int k = 0; k < N; k++)
			for (int l = 0; l < M; l++)
			{
				double r = Distance.distance(k-N/2, l-M/2);
				double phi = atan(k-N/2, l-M/2) - lattPhi;
				radial[k][l] = r*Math.exp(-r/L);
				smooth[k][l] = radial[k][l]*Math.sin(4*phi);
			}
		double[][] smoothAns = new double[N][M];
		FieldOps.shift(smooth, smoothAns);
		
		double sum = FieldOps.sum(radial);		
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				smoothAns[i][j]/=sum;
		
		return smoothAns;
	}
	public static double[][] gaussSmooth(double[][] data, double li, double lj)
	{
		int gli = (int)(3*(li+1));
		int glj = (int)(3*(lj+1));
		int N = data.length;
		int M = data[0].length;
		double[][] smooth = new double[N][M];
		
		if (li == 0 || lj == 0){
			return data.clone();
		}
		if (li == Double.POSITIVE_INFINITY || li == -1)
		{
			double d = mean(data);
			for (int i = 0; i < N; i++)//if (i % 10 == 0) System.out.println(i + " out of " + size);
				for (int j = 0; j < M; j++)
					smooth[i][j] = d;
			return smooth;
		}
		
		double[][] gauss = new double [gli+1][glj+1];

		for (int k = 0; k < gli+1; k++)
			for (int l = 0; l < glj+1; l++)
				gauss[k][l] = Math.exp(-((k*k)/(2*li*li) + (l*l)/(2*lj*lj)));
		
		double gsum;
		int xprime, yprime, xpmin, xpmax, ypmin, ypmax;
		for (int i = 0; i < N; i++){//if (i % 10 == 0) System.out.println(i + " out of " + size);
			for (int j = 0; j < M; j++)
			{
				gsum = 0;
				smooth[i][j] = 0;
				//this is done so that xprime runs from 0 to N, or from x - biggerL/2 to x+biggerL/2
				xpmin = Math.max(0, i - gli);
				xpmax = Math.min(N, i + gli);
				ypmin = Math.max(0, j - glj);
				ypmax = Math.min(M, j + glj);
				for (xprime = xpmin; xprime < xpmax; xprime++)
					for (yprime = ypmin; yprime < ypmax; yprime++)
					{
						gsum += gauss[Math.abs(i-xprime)][Math.abs(j-yprime)];
						smooth[i][j] += gauss[Math.abs(i-xprime)][Math.abs(j-yprime)]*data[xprime][yprime];
					}
				smooth[i][j] /= gsum;
			}
		}
		return smooth;
	}
	public static double[][] subtractGaussSmooth(double[][] data, double L)
	{
		double[][] gausSmooth = gaussSmooth(data, L);
		return FieldOps.minus(data, gausSmooth);
	}
	public static double[][][] subtractGaussSmooth_ReturnBoth(double[][] data, double L)
	{
		double[][] gausSmooth = gaussSmooth(data, L);
		return new double[][][] {FieldOps.minus(data, gausSmooth), gausSmooth};
	}
	public static double[][] cutOffExtremes(double[][] data, double[][] localAvg, double gamma)
	{
		double sigma = FieldOps.sigma(data);
		double[][] answer = new double [data.length][data[0].length];
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
			{
				answer[i][j] = Math.max(data[i][j], localAvg[i][j] - gamma*sigma);
				answer[i][j] = Math.min(answer[i][j], localAvg[i][j] + gamma*sigma);
			}
		return answer;
	}
	public static double[][] cutOffExtremes(double[][] data, double[][] localAvg, double gammaup, double gammadown, double sigma)
	{
//		sigma = FieldOps.sigma(data);
		double[][] answer = new double [data.length][data[0].length];
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
			{
				answer[i][j] = Math.max(data[i][j], localAvg[i][j] - gammadown*sigma);
				answer[i][j] = Math.min(answer[i][j], localAvg[i][j] + gammaup*sigma);
			}
		return answer;
	}
	public static double[][] cutOffExtremes(double[][] data, double gammaup, double gammadown, double sigma)
	{
//		sigma = FieldOps.sigma(data);
		double[][] answer = new double [data.length][data[0].length];
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
			{
				answer[i][j] = Math.max(data[i][j], 0 - gammadown*sigma);
				answer[i][j] = Math.min(answer[i][j], 0 + gammaup*sigma);
			}
		return answer;
	}
	public static void cutOffExtremes(double[][] data, double[][] localAvg, double gammaup, double gammadown, double sigma, double[][] answer)
	{
		//gammadown may now be negative.
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
			{
				answer[i][j] = Math.max(data[i][j], localAvg[i][j] + gammadown*sigma);
				answer[i][j] = Math.min(answer[i][j], localAvg[i][j] + gammaup*sigma);
			}
	}
	public static void cutOffExtremes(double[][] data, double[][] localAvg, double gammaup, double gammadown, double sigma, double[][] answer, int[][] record)
	{
		//gammadown may now be negative.
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
			{
				record[i][j] = 0;
				answer[i][j] = Math.max(data[i][j], localAvg[i][j] + gammadown*sigma);
				answer[i][j] = Math.min(answer[i][j], localAvg[i][j] + gammaup*sigma);

				if (answer[i][j] < data[i][j]) record[i][j] = 1;
				else if (answer[i][j] > data[i][j]) record[i][j] = -1;
			}
	}
	public static double[][] replaceExtremes(double[][] data, double[][] localAvg, double gammaup, double gammadown, double sigma)
	{
		double[][] answer = new double [data.length][data[0].length];
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
				if (data[i][j] > localAvg[i][j] + gammaup*sigma || data[i][j] < localAvg[i][j] - gammadown*sigma)
					answer[i][j] = localAvg[i][j];
				else
					answer[i][j] = data[i][j];
		return answer;
	}
	public static double[][] replaceExtremes(double[][] data, double gammaup, double gammadown, double sigma, double gammarep)
	{
		double[][] answer = new double [data.length][data[0].length];
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
				if (data[i][j] > gammaup*sigma || data[i][j] < - gammadown*sigma)
					answer[i][j] = gammarep*sigma;
				else
					answer[i][j] = data[i][j];
		return answer;
	}
	public static void cutOffExtremes(double[][] data, double[][] localAvg, double gamma, double[][] answer)
	{
		double sigma = FieldOps.sigma(data);
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
			{
				answer[i][j] = Math.max(data[i][j], localAvg[i][j] - gamma*sigma);
				answer[i][j] = Math.min(answer[i][j], localAvg[i][j] + gamma*sigma);
			}
	}
	
//	public static void applyLinearTransformation(double[][] data, double[][] output)
	public static double[][][] getU(double[][] phasex, double[][] phasey, double[][] bragg, double eta)
	{
		int N = phasex.length;
		int M = phasex.length;
		int sizediff = N/M;
		if (sizediff == 2)
		{
			phasex = expand(phasex);
			phasey = expand(phasey);
		}
		if (sizediff == 4)
		{
			phasex = expand4(phasex);
			phasey = expand4(phasey);
		}
		double[][][] u = new double [N][N][2];
		double[] mags = new double [2]; mags[0] = Complex.mag(bragg[0]); mags[1] = Complex.mag(bragg[1]);
		double[] k1Hat = new double [2], k2Hat = new double [2];
		k1Hat[0] = bragg[0][0]/mags[0]; k1Hat[1] = bragg[0][1]/mags[0];
		k2Hat[0] = bragg[1][0]/mags[1]; k2Hat[1] = bragg[1][1]/mags[1];

		double costh = k1Hat[0]*k2Hat[0] + k1Hat[1]*k2Hat[1];

		//use the gramm-schmidt procedure to form k2PrimeHat from k1Hat and k2Hat.
		double[] k2PrimeHat = new double [2];
		k2PrimeHat[0] = k2Hat[0] - costh*k1Hat[0];
		k2PrimeHat[1] = k2Hat[1] - costh*k1Hat[1];
		double magPrime = Complex.mag(k2PrimeHat);
		k2PrimeHat[0] /= magPrime; k2PrimeHat[1] /= magPrime;
		
		//The dot products of u with the unit vectors
		double udot1, udot2, udot2Prime;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				udot1 = (phasex[i][j]/mags[0]);
				udot2 = (phasey[i][j]/mags[1]);
				udot2Prime = (udot2 - costh*udot1)/magPrime;
				
				u[i][j][0] = eta*(udot1*k1Hat[0] + udot2Prime*k2PrimeHat[0]);
				u[i][j][1] = eta*(udot1*k1Hat[1] + udot2Prime*k2PrimeHat[1]);
			}
		return u;
	}
	/**
	 * returns two vectors, the first being aHat and the second being the Graham-schmidt orthonormal bHat.
	 * @param a
	 * @param b
	 * @return
	 */
	public static double[][] getGrahamSchmidtUnitVectors(double[] a, double[] b)
	{
		double[] mags = new double [2]; mags[0] = Complex.mag(a); mags[1] = Complex.mag(b);
		double[] k1Hat = new double [2], k2Hat = new double [2];
		k1Hat[0] = a[0]/mags[0]; k1Hat[1] = a[1]/mags[0];
		k2Hat[0] = b[0]/mags[1]; k2Hat[1] = b[1]/mags[1];

		double costh = k1Hat[0]*k2Hat[0] + k1Hat[1]*k2Hat[1];

		//use the gramm-schmidt procedure to form k2PrimeHat from k1Hat and k2Hat.
		double[] k2PrimeHat = new double [2];
		k2PrimeHat[0] = k2Hat[0] - costh*k1Hat[0];
		k2PrimeHat[1] = k2Hat[1] - costh*k1Hat[1];
		double magPrime = Complex.mag(k2PrimeHat);
		k2PrimeHat[0] /= magPrime; k2PrimeHat[1] /= magPrime;
		return new double[][] {k1Hat, k2PrimeHat};

	}
	public static void putU (double[][] phasex, double[][] phasey, double[][] bragg, double eta, double[][][] u)
	{
		int N = phasex.length;
		int M = phasex.length;
		int sizediff = N/M;
		if (sizediff == 2)
		{
			phasex = expand(phasex);
			phasey = expand(phasey);
		}
		if (sizediff == 4)
		{
			phasex = expand4(phasex);
			phasey = expand4(phasey);
		}
		double[] mags = new double [2]; mags[0] = Complex.mag(bragg[0]); mags[1] = Complex.mag(bragg[1]);
		double[] k1Hat = new double [2], k2Hat = new double [2];
		k1Hat[0] = bragg[0][0]/mags[0]; k1Hat[1] = bragg[0][1]/mags[0];
		k2Hat[0] = bragg[1][0]/mags[1]; k2Hat[1] = bragg[1][1]/mags[1];

		double costh = k1Hat[0]*k2Hat[0] + k1Hat[1]*k2Hat[1];

		//use the gramm-schmidt procedure to form k2PrimeHat from k1Hat and k2Hat.
		double[] k2PrimeHat = new double [2];
		k2PrimeHat[0] = k2Hat[0] - costh*k1Hat[0];
		k2PrimeHat[1] = k2Hat[1] - costh*k1Hat[1];
		double magPrime = Complex.mag(k2PrimeHat);
		k2PrimeHat[0] /= magPrime; k2PrimeHat[1] /= magPrime;
		
		//The dot products of u with the unit vectors
		double udot1, udot2, udot2Prime;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				udot1 = (phasex[i][j]/mags[0]);
				udot2 = (phasey[i][j]/mags[1]);
				udot2Prime = (udot2 - costh*udot1)/magPrime;
				
				u[i][j][0] = eta*(udot1*k1Hat[0] + udot2Prime*k2PrimeHat[0]);
				u[i][j][1] = eta*(udot1*k1Hat[1] + udot2Prime*k2PrimeHat[1]);
			}
		}
	
	//in this code we strictly adhere to the custom of regarding the "center" of a pixel as being at (i+0.5, j+0.5) for purposes of bilinear interpolation (last line only)
	//the matrix multiplication procedure was corrected on July 31st 2012.
	public static void applyLinearTransformation(double[][] data, double[][] matrix, double[] origin, double[][] result)
	{
		applyLinearTransformation(data, matrix, origin, result, FieldOps.mean(data));
	}
	
	/**
	 * This returns the integral of v dot dl along the curve contained in the array list, which consists of 
	 * hopefully contiguous pixels in the order of integration. If contour is meant to be a closed curve, the first
	 * and last entry must be the same.
	 * @param vector2nm
	 * @param contour
	 * @return
	 */
	public static double getLineIntegral(double[][][] vector2nm, ArrayList<int[]> contour){
		double sum = 0;
		int[] dr = new int [2];
		//We sum not over the pixels in contour, but on the intervals between them.
		//The value of vector2nm along the interval is taken to be the average value at start and end 
		int x0, x1, y0, y1;
		double[] vector = new double [2];
//		System.out.println(contour.size());
		for (int i = 0; i < contour.size()-1; i++)
		{
			x0 = contour.get(i)[0]; y0 = contour.get(i)[1];
			x1 = contour.get(i+1)[0]; y1 = contour.get(i+1)[1];
			dr[0] = x1-x0;
			dr[1] = y1-y0;
			vector[0] = (vector2nm[0][x0][y0] + vector2nm[0][x1][y1])/2;
			vector[1] = (vector2nm[1][x0][y0] + vector2nm[1][x1][y1])/2;
			sum += Distance.dot(vector[0], vector[1], dr[0], dr[1]);
//			System.out.print(Distance.dot(vector[0], vector[1], dr[0], dr[1]) + "\t");
		}
//		System.out.println();
		return sum;
	}
	
	public static ArrayList<int[]> getSquareContour(int[] center, int halfSide)
	{
		ArrayList<int[]> ans = new ArrayList<int[]>();
		for (int i = -halfSide; i < halfSide; i++)
			ans.add(new int [] {center[0]+halfSide,center[1]+i}); //The +x side, ending at +y.
		for (int i = -halfSide; i < halfSide; i++)
			ans.add(new int [] {center[0]-i,center[1]+halfSide}); //The +y side, ending at -x.
		for (int i = -halfSide; i < halfSide; i++)
			ans.add(new int [] {center[0]-halfSide,center[1]-i}); //The -x side, ending at -y.
		for (int i = -halfSide; i <= halfSide; i++)
			ans.add(new int [] {center[0]+i,center[1]-halfSide}); //The -y side, ending at +x. and repeating te first pixel.
		
		return ans;
	}
	
	public static void applyLinearTransformation(double[][] data, double[][] matrix, double[] origin, double[][] result, double mean)
	{
		int N = data.length;
		int x, y;
		int ox, oy;
		ox = (int)Math.round(origin[0]);
		oy = (int)Math.round(origin[1]);
		double xp, yp;
		double ip, jp;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				result[i][j] = 0;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				x = i-ox;
				y = j-oy;
				xp = matrix[0][0]*x + matrix[1][0]*y;
				yp = matrix[0][1]*x + matrix[1][1]*y;
				
				ip = xp+origin[0];
				jp = yp+origin[1];
				
				result[i][j] = getValueAt(data, ip+0.5, jp+0.5, mean);
				//just set the value at [i][j] equal to the value of the target at the transformed [i'][j'].
				//The latter to be determined by bilinear interpolation since i' and j' are not integers.
				
//				translateOnePixelAlt(data, result, i, j, (xp-x), (yp-y));
			}
	}
	public static double[][][] generateLinearTransformation(double[][] data, double[][] matrix, int[] origin)
	{
		int N = data.length;
		double[][][] trans = new double [N][N][2];
		int x, y;
		double xp, yp;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				x = i-origin[0];
				y = j-origin[1];
				xp = matrix[0][0]*x + matrix[1][0]*y;
				yp = matrix[1][0]*x + matrix[1][1]*y;
				trans[i][j][0] = xp-x; trans[i][j][1] = yp-y;
			}
		return trans;
	}
	
	public static void zero(double[][] x)
	{
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
				x[i][j] = 0;
	}
	public static void zero(double[][] x, boolean[][] doIt)
	{
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
				if (doIt[i][j])
					x[i][j] = 0;
	}
	public static void setEqualTo(double[][] x, double v)
	{
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
				x[i][j] = v;
	}
	public static void zero(double[][] x, int imin, int imax, int jmin, int jmax)
	{
		for (int i = imin; i < imax; i++)
			for (int j = jmin; j < jmax; j++)
				x[i][j] = 0;
	}
	
	public static void autocorrelate(double[][] source, double[][] target)
	{
		int N = source.length, M = source[0].length;
		int xp, yp;
		for (int i = 0; i < N; i++){System.out.print(" " + i);
			for (int j = 0; j < M; j++)
			{
				for (int k = -N; k < 2*N; k++)
					for (int m = -M; m < 2*M; m++)
//				for (int k = 0; k < N; k++)
//					for (int m = 0; m < M; m++)
					{
						xp = i+k;
						yp = j+m;
						if (xp < N && yp < M && xp >= 0 && yp >= 0 && k > 0 && k < M && m > 0 && m < M)
							target[i][j] += source[xp][yp]*source[k][m];
					}
			}
		}
	}
	
	
	//The next couple of methods attempt to calculate the joint density of states
	public static double shiftedSelfSum(double[][] source, int k, int m)
	{
		double ans = 0;
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source.length; j++)
			{
				if (withinBounds(i+k, j+m, source))
					ans += source[i][j]*source[i+k][j+m];
			}
		return ans;
	}
	public static double[][] getJDOS_sameSize(double[][] source)
	{
		int shiftx, shifty;
		double[][] ans = new double [source.length][source[0].length];
		double[][] copy = FieldOps.copy(source);
		FieldOps.minus(copy, FieldOps.mean(source));
		for (int i = 0; i < source.length; i++){System.out.print(" " + i);
			for (int j = 0; j < source.length; j++)
			{
				shiftx =  i-source.length/2;
				shifty = j-source.length/2;
				ans[i][j] = shiftedSelfSum(copy, shiftx, shifty);
			}
		}
		return ans;
	}
	/**
	 * The logic of this method is from a Wikipedia article.
	 * @param source
	 * @return
	 */
	public static double[][] getAutocorrelationFourier(double[][] source)
	{
		int nx = source.length, ny = source[0].length;
		double[][][] ans = new double [source.length][source[0].length][2];
//		double[][] copy = FieldOps.copy(source);
//		FieldOps.minus(copy, FieldOps.mean(source));
		
		double mean = FieldOps.mean(source);
//		FieldOps.minusEquals(source, mean);
		
		double[][][] fftz = new double [nx][ny][2];
		FFTOps.putFFT(source, fftz, false);
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				fftz[i][j][0] = fftz[i][j][0]*fftz[i][j][0] + fftz[i][j][1]*fftz[i][j][1];
				fftz[i][j][1] = 0;
			}
		FFTOps.putIFFT(fftz, ans, false);
		double[][] ret = FieldOps.getIndex(ans, 0);
		ans = null;
		double[][] copy = new double [nx][ny];
		FieldOps.shift(ret, copy, nx/2, ny/2);
//		FieldOps.plusEquals(source, mean);
		return copy;
	}
	public static double[][] getCrosscorrelationFourier(double[][] source, double[][] source2)
	{
		int nx = source.length, ny = source[0].length;
		double[][][] ans = new double [source.length][source[0].length][2];
//		double[][] copy = FieldOps.copy(source);
//		FieldOps.minus(copy, FieldOps.mean(source));
		
		double mean = FieldOps.mean(source);
		FieldOps.minusEquals(source, mean);
		double mean2 = FieldOps.mean(source2);
		FieldOps.minusEquals(source2, mean2);
		
		double[][][] fftz = new double [nx][ny][2];
		double[][][] fftz2 = new double [nx][nx][2];
		double[][][] fftzans = new double [nx][ny][2];
		
		FFTOps.putFFT(source, fftz, false);
		FFTOps.putFFT(source2, fftz2, false);
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				fftzans[i][j][0] = fftz[i][j][0]*fftz2[i][j][0] + fftz[i][j][1]*fftz2[i][j][1];
				fftzans[i][j][1] = -fftz[i][j][0]*fftz2[i][j][1] + fftz2[i][j][0]*fftz[i][j][1];
			}
		FFTOps.putIFFT(fftzans, ans, false);
		double[][] ret = FieldOps.getIndex(ans, 0);
		ans = null;
//		double[][] copy = new double [nx][ny];
//		FieldOps.shift(ret, copy, nx/2, ny/2);
		FieldOps.plusEquals(source, mean);
		FieldOps.plusEquals(source2, mean2);
		return ret;
	}
	
	//These two assume that both data sets are the same size and the sum is over all points
	public static double correlation(double[][] a, double[][] b)
	{
		double[][] a0 = FieldOps.copy(a), b0 = FieldOps.copy(b);
		FieldOps.subtractAvg(a0); FieldOps.subtractAvg(b0);
		return correlationMeanZero(a0, b0);
	}
	public static double correlationMeanZero(double[][] a, double[][] b)
	{
		double[] sigma = new double [2];
		sigma[0] = sigma(a, 0);
		sigma[1] = sigma(b, 0);
		
		double r = 0;
		for (int i = 0; i < a.length; i++)
			for (int j = 0; j < a[0].length; j++)
				r += a[i][j]*b[i][j];
		
		return r/(sigma[0]*sigma[1]*a.length*a[0].length);
	}
	public static double innerProduct(double[][] a, double[][] b)
	{
		double r = 0;
		for (int i = 0; i < a.length; i++)
			for (int j = 0; j < a[0].length; j++)
				r += a[i][j]*b[i][j];
		
		return r;
	}
	
	//Here b is larger in size than a. a is correlated with the region of b starting at [boffi, boffj].
	//It is assumed that the average of a and the LOCAL average of b are zero
	public static double correlationMeanZero(double[][] a, double[][] b, int boffi, int boffj)
	{
		double[] sigma = new double [2];
		sigma[0] = sigma(a, 0);
		sigma[1] = sigmaLocal(b, 0, boffi, boffj, a.length, a[0].length);
		
		double r = 0;
		int ninsum = 0;
		for (int i = 0; i < a.length; i++)
			for (int j = 0; j < a[0].length; j++)
				if (i+boffi < b.length && j+boffj < b[0].length && i+boffi >= 0 && j+boffj >= 0)
				{
					r += a[i][j]*b[i+boffi][j+boffj];
					ninsum++;
				}
				//else the value of b there is assumed to be zero
		if (ninsum > 0)
			return r/(sigma[0]*sigma[1]*ninsum);
		else
			return 0;
	}
	public static double correlation(double[][] a, double[][] b, int boffi, int boffj)
	{
		double[] sigma = new double [2];
		double[][] a0 = FieldOps.copy(a);
		FieldOps.subtractAvg(a0);
		sigma[0] = sigma(a0, 0);
//		sigma[1] = sigmaLocal(b, 0, boffi, boffj, a.length, a[0].length);
		int[] bounds = new int [4];
		bounds(boffi, boffi+a.length, boffj, boffj+a[0].length, b.length, b[0].length, bounds);
		double[][] bsub = FieldOps.subset(b, bounds[0], bounds[1], bounds[2], bounds[3]);
		sigma[1] = sigma(bsub);
		double bsubavg = FieldOps.mean(bsub);
		FieldOps.minusEquals(b, bsubavg);
		
		double r = 0;
		int ninsum = 0;
		for (int i = 0; i < a.length; i++)
			for (int j = 0; j < a[0].length; j++)
				if (i+boffi < b.length && j+boffj < b[0].length && i+boffi >= 0 && j+boffj >= 0)
				{
					r += a[i][j]*b[i+boffi][j+boffj];
					ninsum++;
				}
		
		FieldOps.plusEquals(b, bsubavg);
				//else the value of b there is assumed to be zero
		if (ninsum > 0)
			return r/(sigma[0]*sigma[1]*ninsum);
		else
			return 0;
	}

	public static double getIntegralBetween(double[][] a, double[][] b, int boffi, int boffj)
	{
		double r = 0;
		for (int i = 0; i < a.length; i++)
			for (int j = 0; j < a[0].length; j++)
				if (i+boffi < b.length && j+boffj < b[0].length && i+boffi >= 0 && j+boffj >= 0)
				{
					r += a[i][j]*b[i+boffi][j+boffj];
				}
		return r;
	}
	
	//subtracts integral*a from t, which is supposed to be a copy of b. The integral 
	public static void subtractBAfterIntegral(double[][] a, double[][] b, int boffi, int boffj, double[][] t)
	{
		double integral = getIntegralBetween(a, b, boffi, boffj);
		for (int i = 0; i < a.length; i++)
			for (int j = 0; j < a[0].length; j++)
				if (i+boffi < b.length && j+boffj < b[0].length && i+boffi >= 0 && j+boffj >= 0)
				{
					t[i+boffi][j+boffj] -= integral*a[i][j];
				}

	}
	//a is smaller than b;
	public static double[][] getCorrelationMap(double[][] a, double[][] b)
	{
		double[][] answer = new double [b.length-a.length][b[0].length-a[0].length];
		double[][] a0 = FieldOps.copy(a), b0 = FieldOps.copy(b);
//		FieldOps.subtractAvg(a0); FieldOps.subtractAvg(b0);
		FieldOps.subtractParabolicFit(a0); FieldOps.subtractParabolicFit(b0);
		System.out.print("Done (out of " + answer.length + "): ");
		for (int i = 0; i < answer.length; i++){
			for (int j = 0; j < answer[0].length; j++)
			{
				answer[i][j] = correlationMeanZero(a0, b0, i, j);
			}
			System.out.print("" + (i+1) + " ");
		}
		System.out.println();
		return answer;
	}
	public static double[] getCorrelationvsTheta(double[][] a, double[][] b, int npts, double thmax)
	{
		double[] answer = new double [npts];
		double[][] a0 = FieldOps.copy(a), b0 = FieldOps.copy(b);
		FieldOps.subtractAvg(a0); FieldOps.subtractAvg(b0);
		double[][] target = FieldOps.copy(b0);
		double[][] matrix = new double [2][2];
		double theta = 0;
		double[] origin = new double [] {a.length/2, a[0].length/2};
		double dt = thmax/npts;
		//		System.out.print("Done (out of " + answer.length + "): ");
		for (int i = 0; i < npts; i++){
				theta = i*dt;
				Matrix.putRotationMatrix(theta, matrix);
				FieldOps.applyLinearTransformation(b0, matrix, origin, target);
				answer[i] = correlationMeanZero(a0, target);
			}
		return answer;
	}
	public static void putCorrelationvsTheta(double[][] a, double[][] b, double thmax, double[] answer)
	{
		int npts = answer.length;
		double[][] a0 = FieldOps.copy(a), b0 = FieldOps.copy(b);
		FieldOps.subtractAvg(a0); FieldOps.subtractAvg(b0);
		double[][] target = FieldOps.copy(b0);
		double[][] matrix = new double [2][2];
		double theta = 0;
		double[] origin = new double [] {a.length/2, a[0].length/2};
		double dt = thmax/npts;
		//		System.out.print("Done (out of " + answer.length + "): ");
		for (int i = 0; i < npts; i++){
				theta = i*dt;
				Matrix.putRotationMatrix(theta, matrix);
				FieldOps.applyLinearTransformation(b0, matrix, origin, target);
				answer[i] = correlationMeanZero(a0, target);
			}
	}
	
	/**
	 * The vorticity is defined as follows. Take a path around the pixel (8 nearest neighbors). If from one pixel to the next goes up by more than the threshold
	 * value, add +1. If it goes down, add -1. This way the +1 and -1 will cancel along lines of 2pi change, but not at the points where these lines end.
	 * 
	 * The "line" will end at a corner shared by four pixels. All four pixels will have the +1 or -1 in them. (Unless they are at the edge).
	 * @param data
	 * @param target
	 */
	public static void putVorticity(double[][] data, double[][] target, double threshold)
	{
		for (int i = 0; i < target.length; i++)
			for (int j = 0; j < target[0].length; j++)
			{
				target[i][j] = 0;
				if (i > 0 && i < target.length-1 && j > 0 && j < target[0].length-1)
				{
					if (data[i+1][j] - data[i+1][j+1] > threshold) target[i][j]++;
					else if (data[i+1][j] - data[i+1][j+1] < -threshold) target[i][j]--;
					if (data[i+1][j+1] - data[i][j+1] > threshold) target[i][j]++;
					else if (data[i+1][j+1] - data[i][j+1] < -threshold) target[i][j]--;
					if (data[i][j+1] - data[i-1][j+1] > threshold) target[i][j]++;
					else if (data[i][j+1] - data[i-1][j+1] < -threshold) target[i][j]--;
					if (data[i-1][j+1] - data[i-1][j] > threshold) target[i][j]++;
					else if (data[i-1][j+1] - data[i-1][j] < -threshold) target[i][j]--;
					if (data[i-1][j] - data[i-1][j-1] > threshold) target[i][j]++;
					else if (data[i-1][j] - data[i-1][j-1] < -threshold) target[i][j]--;
					if (data[i-1][j-1] - data[i][j-1] > threshold) target[i][j]++;
					else if (data[i-1][j-1] - data[i][j-1] < -threshold) target[i][j]--;
					if (data[i][j-1] - data[i+1][j-1] > threshold) target[i][j]++;
					else if (data[i][j-1] - data[i+1][j-1] < -threshold) target[i][j]--;
					if (data[i+1][j-1] - data[i+1][j] > threshold) target[i][j]++;
					else if (data[i+1][j-1] - data[i+1][j] < -threshold) target[i][j]--;
				}
			}
	}
	
	public static double[][] getVorticity(double[][] source, double threshold)
	{
		double[][] ans = new double[source.length][source[0].length];
		putVorticity(source, ans, threshold);
		return ans;
	}
	public static double getNVortices(double[][] data, double threshold)
	{
		double[][] vort = getVorticity(data, threshold);
		FieldOps.abs(vort);
		return ArrayOps.sum(vort)/4;
	}
	
	/**
	 * returns the number of times adjacent pixels differ from each other by more than "Threshold"
	 * @param data
	 * @param threshold
	 * @return
	 */
	public static double getNLargeDifferences(double[][] f, double threshold)
	{
		double ans = 0;
		int nx = f.length, ny = f[0].length;
		//the four corner pixels each have two neighbors:
		if (Math.abs(f[0][0] - f[0][1]) > threshold)
			ans++;
		if (Math.abs(f[0][0] - f[1][0]) > threshold)
			ans++;
		if (Math.abs(f[nx-1][0] - f[nx-1][1]) > threshold)
			ans++;
		if (Math.abs(f[nx-1][0] - f[nx-2][0]) > threshold)
			ans++;
		if (Math.abs(f[0][ny-1] - f[0][ny-2]) > threshold)
			ans++;
		if (Math.abs(f[0][ny-1] - f[1][ny-1]) > threshold)
			ans++;
		if (Math.abs(f[nx-1][ny-1] - f[nx-1][ny-2]) > threshold)
			ans++;
		if (Math.abs(f[nx-1][ny-1] - f[nx-2][ny-1]) > threshold)
			ans++;
		
		int i = 0, j = 0;
		//the four edges each have 3 neighbors:
		for (i = 1; i < ny-1; i++)
		{
			if (Math.abs(f[1][i] - f[0][i]) > threshold) ans++;
			if (Math.abs(f[nx-1][i] - f[nx-2][i]) > threshold) ans++;
		}
		for (i = 1; i < nx-1; i++)
		{
			if (Math.abs(f[i][1] - f[i][0]) > threshold) ans++;
			if (Math.abs(f[i][ny-1] - f[i][ny-2]) > threshold) ans++;
		}
		
		for (i = 1; i < nx - 1; i++)
			for (j = 1; j < ny - 1; j++)
			{
				if (Math.abs(f[i][j] - f[i+1][j]) > threshold) ans++;
				if (Math.abs(f[i][j] - f[i-1][j]) > threshold) ans++;
				if (Math.abs(f[i][j] - f[i][j+1]) > threshold) ans++;
				if (Math.abs(f[i][j] - f[i][j-1]) > threshold) ans++;
			}

		return ans;
	}
	//a is smaller than b; returns an object the size of b
	public static double[][] getCorrelationMapCent(double[][] a, double[][] b)
	{
		double[][] answer = new double [b.length][b[0].length];
		double[][] a0 = FieldOps.copy(a), b0 = FieldOps.copy(b);
		FieldOps.subtractAvg(a0); FieldOps.subtractAvg(b0);
		double mean;
		//subtract local average from a0 and re-add it.
//		System.out.print("Done (out of " + answer.length + "): ");
		for (int i = 0; i < answer.length; i++){
			for (int j = 0; j < answer[0].length; j++)
			{
				mean = FieldOps.subtractAvg(b0, i-a.length/2, j-a[0].length/2, a.length, a[0].length);
				answer[i][j] = correlationMeanZero(a0, b0, i-a.length/2, j-a[0].length/2);
				FieldOps.addValue(b0, i-a.length/2, j-a[0].length/2, a.length, a[0].length, mean);
			}
//			if (i % 50 == 0) System.out.println();
//			System.out.print("" + (i+1) + " ");
		}
//		System.out.println();
		return answer;
	}
	public static double[][] getIntegralMapCent(double[][] a, double[][] b)
	{
		double[][] answer = new double [b.length][b[0].length];
		double[][] a0 = FieldOps.copy(a), b0 = FieldOps.copy(b);
		FieldOps.subtractAvg(a0); FieldOps.subtractAvg(b0);
		double mean;
		//subtract local average from a0 and re-add it.
//		System.out.print("Done (out of " + answer.length + "): ");
		for (int i = 0; i < answer.length; i++){
			for (int j = 0; j < answer[0].length; j++)
			{
				answer[i][j] = FieldOps.getIntegralBetween(a, b, i-a.length/2, j-a[0].length/2);
			}
//			if (i % 50 == 0) System.out.println();
//			System.out.print("" + (i+1) + " ");
		}
//		System.out.println();
		return answer;
	}
	
	public static double[] dump(double[][] source)
	{
		double[] dump = new double [source.length*source[0].length];
		int n = 0;
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
				dump[n++] = source[i][j];
		return dump;
	}
	public static double[] dump(double[][][] source)
	{
		double[] dump = new double [source[0][0].length*source.length*source[0].length];
		int n = 0;
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[i].length; j++)
				for (int k = 0; k < source[i][j].length; k++)
					dump[n++] = source[i][j][k];
		return dump;
	}
	public static double[] getSortedDump(double[][] source)
	{
		double[] dump = dump(source);
		ArrayOps.quicksort(dump);
		return dump;
	}
	public static double[] getSortedDump(double[][][] source)
	{
		double[] dump = dump(source);
		ArrayOps.quicksort(dump);
		return dump;
	}
	public static boolean[][] isGreaterThan(double[][] source, double cutoff)
	{
		boolean[][] ans = new boolean [source.length][source[0].length];
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
				ans[i][j] = source[i][j] > cutoff;
		return ans;
	}
	public static boolean[][] getAnd(boolean[][] a, boolean[][] b)
	{
		boolean[][] ans = new boolean [a.length][a[0].length];
		for (int i = 0; i < a.length; i++)
			for (int j = 0; j < a[0].length; j++)
				ans[i][j] = a[i][j] && b[i][j];
		return ans;
	}
	public static void isGreaterThan(double[][] source, double cutoff, boolean[][] ans)
	{
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
				ans[i][j] = source[i][j] > cutoff;
	}
	public static void cutOffExtremes(double[][] source, double lowerperc, double upperperc)
	{
		double[] sorted = getSortedDump(source);
		int n = sorted.length;
		double min = sorted[(int)(lowerperc*(n-1))];
		double max = sorted[(int)(upperperc*(n-1))];
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
			{
				source[i][j] = Math.max(source[i][j], min);
				source[i][j] = Math.min(source[i][j], max);
			}
	}
	public static double[][][] cutOffExtremes3D(double[][][] source, double lowerperc, double upperperc, boolean kill, boolean returnAnything)
	{
		double[] sorted = getSortedDump(source);
		int n = sorted.length;
		double min = sorted[(int)(lowerperc*(n-1))];
		double max = sorted[(int)(upperperc*(n-1))];
		double[][][] oldSource = null;
		if (returnAnything) oldSource = FieldOps.copy(source);
		if (!kill)
			for (int i = 0; i < source.length; i++)
				for (int j = 0; j < source[i].length; j++)
					for (int k = 0; k < source[i][j].length; k++)
					{
						source[i][j][k] = Math.max(source[i][j][k], min);
						source[i][j][k] = Math.min(source[i][j][k], max);
					}
		else
			for (int i = 0; i < source.length; i++)
				for (int j = 0; j < source[i].length; j++)
					for (int k = 0; k < source[i][j].length; k++)
						if (source[i][j][k] > max || source[i][j][k] < min) source[i][j][k] = 0; else;
	
		if (returnAnything) FieldOps.minusEquals(oldSource, source);
		return oldSource;
	
	}
	/**
	 * This cuts off the extremes of the source (modifies source), and returns a double array of -1, 0, or 1 depending if the pixel was too high, OK, or too low. 
	 * @param source
	 * @param lowerperc
	 * @param upperperc
	 * @return
	 */
	public static double[][] cutOffExtremesPixels(double[][] source, double lowerperc, double upperperc)
	{
		double[] sorted = getSortedDump(source);
		double old;
		double[][] ans = new double [source.length][source[0].length];
		int n = sorted.length;
		double min = sorted[(int)(lowerperc*n)];
		double max = sorted[(int)(upperperc*n)];
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
			{
				old = source[i][j];
				source[i][j] = Math.max(source[i][j], min);
				source[i][j] = Math.min(source[i][j], max);
				if (source[i][j] < old)
					ans[i][j]  = 1;
				else if (source[i][j] > old)
					ans[i][j] = -1;
			}
		return ans;
	}
	public static void replaceExtremes(double[][] source, double lowerperc, double upperperc, double value)
	{
		double[] sorted = getSortedDump(source);
		int n = sorted.length;
		double min = sorted[(int)(lowerperc*n)];
		double max = sorted[(int)(upperperc*n)];
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
			{
				if (source[i][j] > max || source[i][j] < min) source[i][j] = value;
			}
	}
	public static void cutOffExtremes(double[][] source, double lowerperc, double upperperc, double[] sorted)
	{
		int n = sorted.length;
		double min = sorted[(int)(lowerperc*n)];
		double max = sorted[(int)(upperperc*n)];
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
			{
				source[i][j] = Math.max(source[i][j], min);
				source[i][j] = Math.min(source[i][j], max);
			}
	}
	public static void replaceExtremes(double[][] source, double lowerperc, double upperperc, double value, double[] sorted)
	{
		int n = sorted.length;
		double min = sorted[(int)(lowerperc*n)];
		double max = sorted[(int)(upperperc*n)];
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
			{
				if (source[i][j] > max || source[i][j] < min) source[i][j] = value;
			}
	}
	public static void cutOffExtremes(double[][] source, double lowerperc, double upperperc, double[] sorted, double[][] target, int[][] record)
	{
		if (lowerperc < 0 || lowerperc >= 1 || upperperc < 0 || upperperc >= 1) return;
//		System.out.println(lowerperc + ", " + upperperc);
		int n = sorted.length;
		double min = sorted[(int)(lowerperc*n)];
		double max = sorted[(int)(upperperc*n)];
		System.out.println(min + ", " + max);
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
			{
//				target[i][j] = 0;
				record[i][j] = 0;
				target[i][j] = Math.max(source[i][j], min);
				target[i][j] = Math.min(target[i][j], max);
				
				if (target[i][j] < source[i][j]) record[i][j] = 1;
				else if (target[i][j] > source[i][j]) record[i][j] = -1;

			}
	}
	public static void replaceExtremes(double[][] source, double lowerperc, double upperperc, double value, double[] sorted, double[][] target)
	{
		int n = sorted.length;
		double min = sorted[(int)(lowerperc*n)];
		double max = sorted[(int)(upperperc*n)];
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
			{
				if (source[i][j] > max || source[i][j] < min) target[i][j] = value;
			}
	}
	public static void replaceExtremesValue(double[][] source, double min, double max, double value, double[][] target)
	{
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
			{
				if (source[i][j] > max || source[i][j] < min) target[i][j] = value;
				else target[i][j] = source[i][j];
			}
	}
	public static void cutoffExtremesValue(double[][] source, double min, double max, double[][] target)
	{
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
			{
				if (source[i][j] > max) target[i][j] = max;
				else if (source[i][j] < min) target[i][j] = min;
				else target[i][j] = source[i][j];
			}
	}
	public static void cutoffExtremesValue(double[][] source, double min, double max)
	{
		for (int i = 0; i < source.length; i++)
			for (int j = 0; j < source[0].length; j++)
			{
				if (source[i][j] > max) source[i][j] = max;
				else if (source[i][j] < min) source[i][j] = min;
			}
	}
	
	public static double getValueAt(double[][] field, double x, double y, double mean)	{
		//This is with bilinear interpolation. The definition of the center is n+1/2.
		//then the important variables are the distances of the point to the points n-1/2 and n+1/2.
		if (x < 0 || y < 0 || x > field.length || y > field[0].length) return mean;
		
		int xx = (int)x, yy = (int)y;
		int xm = xx, xp = xx, ym = yy, yp = yy;	//initial values are arbitrary.
		double a = x-xx, b = y-yy;
		double ap, bp;
		if (a <= 0.5 && b <= 0.5) {xm = xx-1; xp = xx; ym = yy-1; yp = yy;}
		if (a > 0.5 && b <= 0.5) {xm = xx; xp = xx+1; ym = yy-1; yp = yy;}
		if (a <= 0.5 && b > 0.5) {xm = xx-1; xp = xx; ym = yy; yp = yy+1;}
		if (a > 0.5 && b > 0.5) {xm = xx; xp = xx+1; ym = yy; yp = yy+1;}
		ap = x-(xm+0.5);
		bp = y-(ym+0.5); 
		
		if (xm < 0) {xm = 0; xx = 0;}
		if (ym < 0) {ym = 0; yy = 0;}
		if (xp > field.length-1) {xp = field.length-1; xx = field.length-1;}
		if (yp > field[0].length-1) {yp = field[0].length-1; yy = field[0].length-1;}
		return field[xm][ym]*(1-ap)*(1-bp) + field[xp][ym]*ap*(1-bp) + field[xm][yp]*(1-ap)*bp + field[xp][yp]*ap*bp;
	}
	public static double getValueAt(double[][][] field, double x, double y, int thirdindex)	{
		//This is with bilinear interpolation. The definition of the center is n+1/2.
		int xx = (int)x, yy = (int)y;
		int xm = xx, xp = xx, ym = yy, yp = yy;	//initial values are arbitrary.
		double a = x-xx, b = y-yy;
		double ap = 0, bp = 0;
		if (a <= 0.5 && b <= 0.5) {xm = xx-1; xp = xx; ym = yy-1; yp = yy;}
		if (a > 0.5 && b <= 0.5) {xm = xx; xp = xx+1; ym = yy-1; yp = yy;}
		if (a <= 0.5 && b > 0.5) {xm = xx-1; xp = xx; ym = yy; yp = yy+1;}
		if (a > 0.5 && b > 0.5) {xm = xx; xp = xx+1; ym = yy; yp = yy+1;}
		ap = x-(xm+0.5);
		bp = y-(ym+0.5); 
	
		if (xm < 0) {xm = 0; xx = 0;}
		if (ym < 0) {ym = 0; yy = 0;}
		if (xp > field.length-1) {xp = field.length-1; xx = field.length-1;}
		if (yp > field[0].length-1) {yp = field[0].length-1; yy = field[0].length-1;}
		return field[xm][ym][thirdindex]*(1-ap)*(1-bp) + field[xp][ym][thirdindex]*ap*(1-bp) + field[xm][yp][thirdindex]*(1-ap)*bp + field[xp][yp][thirdindex]*ap*bp;
	}
	
	public static double[][] expandBi(double[][] data, int factor, int xmin, int dx, int ymin, int dy)
	{
		double mean = FieldOps.mean(data, xmin, xmin+dx, ymin, ymin+dy);
		double[] r = new double [2];
		double[][] target = new double [dx*factor][dy*factor];
		double dr = 1/(double)(factor), edge = 1/(double)(2*factor);
		for (int i = 0; i < dx*factor; i++)
			for (int j = 0; j < dy*factor; j++)
			{
				r[0] = xmin + edge + i*dr;
				r[1] = ymin + edge + j*dr;
				target[i][j] = getValueAt(data, r[0], r[1], mean);
			}
		return target;
	}
	public static double[][] expandBi(double[][] data, int factor)
	{
		return expandBi(data, factor, 0, data.length, 0, data[0].length);
	}
	public static void expandBi(double[][] data, int factor, int xmin, int dx, int ymin, int dy, double[][] target)
	{
		double mean = FieldOps.mean(data, xmin, xmin+dx, ymin, ymin+dy);
		double[] r = new double [2];
		double dr = 1/(double)(factor), edge = 1/(double)(2*factor);
		for (int i = 0; i < dx*factor; i++)
			for (int j = 0; j < dy*factor; j++)
			{
				r[0] = xmin + edge + i*dr;
				r[1] = ymin + edge + j*dr;
				target[i][j] = getValueAt(data, r[0], r[1], mean);
			}
	}
	public static void expandBiRotatedAboutCenter(double[][] data, double theta, int factor, int xmin, int dx, int ymin, int dy, double[][] target)
	{
		double[][] matrix = new double [2][2];
		Matrix.putRotationMatrix(theta, matrix);
		double[] cent = new double[] {xmin + dx/2.0, ymin + dy/2.0};
		
		
		double mean = FieldOps.mean(data, xmin, xmin+dx, ymin, ymin+dy);
		double[] r = new double [2];
		double dr = 1/(double)(factor), edge = 1/(double)(2*factor);
		for (int i = 0; i < dx*factor; i++)
			for (int j = 0; j < dy*factor; j++)
			{
				r[0] = xmin + edge + i*dr;
				r[1] = ymin + edge + j*dr;
				Matrix.putProductWithOrigin(matrix, r, r, cent);
				target[i][j] = getValueAt(data, r[0], r[1], mean);
			}
	}
	public static void expandBiRotated(double[][] data, double theta, double[] rotOrigin, int factor, int xmin, int dx, int ymin, int dy, double[][] target)
	{
		double[][] matrix = new double [2][2];
		Matrix.putRotationMatrix(theta, matrix);
	
		double mean = FieldOps.mean(data, xmin, xmin+dx, ymin, ymin+dy);
		double[] r = new double [2];
		double dr = 1/(double)(factor), edge = 1/(double)(2*factor);
		for (int i = 0; i < dx*factor; i++)
			for (int j = 0; j < dy*factor; j++)
			{
				r[0] = xmin + edge + i*dr;
				r[1] = ymin + edge + j*dr;
				Matrix.putProductWithOrigin(matrix, r, r, rotOrigin);
				target[i][j] = getValueAt(data, r[0], r[1], mean);
			}
	}
	
	//if nx = dx and ny = d
	public static void evaluateBi(double[][] data, int nx, int ny, int xmin, int dx, int ymin, int dy, double[][] target)
	{
		double mean = FieldOps.mean(data, xmin, xmin+dx, ymin, ymin+dy);
		double[] r = new double [2];
		double dxe = dx/(double)nx, dye = dy/(double)ny;
		double edgex = dxe/2;
		double edgey = dye/2;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				r[0] = xmin + edgex + i*dxe;
				r[1] = ymin + edgey + j*dye;
				target[i][j] = getValueAt(data, r[0], r[1], mean);
			}
	}
	public static void evaluateBi(double[][] data, int nx, int ny, double[][] target)
	{
		int xmin = 0, ymin = 0;
		int dx = data.length, dy = data[0].length;
		double mean = FieldOps.mean(data, xmin, xmin+dx, ymin, ymin+dy);
		double[] r = new double [2];
		double dxe = dx/(double)nx, dye = dy/(double)ny;
		double edgex = dxe/2;
		double edgey = dye/2;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				r[0] = xmin + edgex + i*dxe;
				r[1] = ymin + edgey + j*dye;
				target[i][j] = getValueAt(data, r[0], r[1], mean);
			}
	}
	public static double evaluateAtInt(double[][] data, double px, double py, double mean)
	{
		if (FieldOps.withinBounds(FieldOps.round(px), FieldOps.round(py), data))
			return data[FieldOps.round(px)][FieldOps.round(py)];
		else return mean;
	}

	public static void expandBi(double[][] data, int factor, double[][] target)
	{
		expandBi(data, factor, 0, data.length, 0, data[0].length, target);
	}
	
	/**
	 * The x and y arrays are in pixel units. The ans array is [x.length][y.length].
	 * The interpolating function is assumed to be not null and to contain everything about the data field.
	 * @param x
	 * @param y
	 * @param data
	 */
	public static void putBicubicInterpolation(double[] x, double[] y, BicubicSplineInterpolatingFunction interp, double[][] ans)
	{
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < y.length; j++)
				ans[i][j] = interp.value(x[i], y[j]);
	}
	public static double[][] convert(boolean[][] f)
	{
		double[][] a = new double [f.length][f[0].length];
		for (int i = 0; i < f.length; i++)
			for (int j = 0; j < f[0].length; j++)
			{
				a[i][j] = f[i][j] ? 1 : 0;
			}
		return a;
	}
	public static void convert(boolean[][] f, double[][] a)
	{
		for (int i = 0; i < f.length; i++)
			for (int j = 0; j < f[0].length; j++)
			{
				a[i][j] = f[i][j] ? 1 : 0;
			}
	}
	
	public static double getLocalAvgCircle(double[][] data, double x, double y, double radius)
	{
		double dx, dy; int npts = 0;
		int x0, x1, y0, y1;
		x0 = (int)Math.max(0, x-radius);
		x1 = (int)Math.min(data.length-1, x+radius);
		y0 = (int)Math.max(0, y-radius);
		y1 = (int)Math.min(data[0].length-1, y+radius);
		
		double ans = 0;
		for (int i = x0; i < x1; i++)
			for (int j = y0; j < y1; j++)
			{
				dx = x-i; dy = y-j;
				if (dx*dx + dy*dy <= radius*radius)
				{
					ans += data[i][j];
					npts++;
				}
			}
		return ans/npts;
	}
	
	public static void bounds(int imin, int imax, int jmin, int jmax, int sizei, int sizej, int[] bounds)
	{
		bounds[0] = Math.max(imin, 0);
		bounds[1] = Math.min(imax, sizei);
		bounds[2] = Math.max(jmin, 0);
		bounds[3] = Math.min(jmax,  sizej);
	}
	
	public static void changeZeroToAverage(double[][] data)
	{
		double mean = ArrayOps.meanExceptZero(data);
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[i].length; j++)
				if (data[i][j] == 0)
					data[i][j] = mean;
	}
	public static void changeValueToAverage(double[][] data, double value)
	{
		double mean = ArrayOps.meanExceptValue(data, value);
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[i].length; j++)
				if (data[i][j] == value)
					data[i][j] = mean;
	}
	public static void changeNaNToAverage(double[][] data)
	{
		String nan = "NaN";
		double mean = ArrayOps.meanExceptValueString(data, nan);
		String s;
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[i].length; j++)
			{
				s = "" + data[i][j];
				if (s.equals(nan))
					data[i][j] = mean;
			}
	}
	public static void changeToAverageAfterIndex(double[][] data, int n)
	{
		double mean = ArrayOps.meanExceptValueString(data, "NaN");
		int ind = 0;
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[i].length; j++)
			{
				if (ind++ > n)
					data[i][j] = mean;
			}
	}
	public static void changeToAverageOutsideBox(double[][] data, int xi, int xf, int yi, int yf)
	{
		double mean = ArrayOps.meanInsideBox(data, xi, xf, yi, yf);
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[i].length; j++)
			{
				if (!(i >= xi && i < xf && j >= yi && j < yf))
					data[i][j] = mean;
			}
	}
	/**
	 * The mask is assumed to go from 0 to 1. A calue of zero means keep the data the same. A value of 1 means, replace the data with the "average".
	 * 
	 * @param data
	 * @param mask
	 */
	public static void changeToWithoutMask(double[][] data, double[][] mask)
	{
		double mean = ArrayOps.meanWithoutMask(data, mask);
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[i].length; j++)
			{
				data[i][j] = mask[i][j]*mean + (1-mask[i][j])*data[i][j];
			}
	}
	public static void putGaussDots(double[][] target, double[][] gaussPts, double spread)
	{
		double g;
		double d;
		FieldOps.zero(target);
		for (int i = 0; i < target.length; i++)
			for (int j = 0; j < target[0].length; j++){
				for (int k = 0; k < gaussPts.length; k++)
				{
					d = Distance.distance(i-gaussPts[k][0], j-gaussPts[k][1]);
					g = Math.exp(-(d*d)/(spread*spread));
					target[i][j] += g;
				}
				target[i][j] = Math.min(target[i][j], 1);
			}
	}
	public static void putGaussDots(double[][] target, double[][][] gaussPts, double spread, boolean[][] include)
	{
		double g;
		double d;
		FieldOps.zero(target);
		for (int i = 0; i < target.length; i++)
			for (int j = 0; j < target[0].length; j++){
				for (int k = 0; k < gaussPts.length; k++)
					for (int m = 0; m < gaussPts[0].length; m++)
						if (include[k][m])
							{
							d = Distance.distance(i-gaussPts[k][m][0], j-gaussPts[k][m][1]);
							if (d < 5*spread){ 
								g = Math.exp(-(d*d)/(spread*spread));
								target[i][j] += g;
							}
				}
				target[i][j] = Math.min(target[i][j], 1);
			}
	}
	public static void putGaussDots(double[][] target, double[][][] gaussPts, double spread, double height, boolean limitTo1)
	{
		double g;
		double d;
		FieldOps.zero(target);
		for (int i = 0; i < target.length; i++){System.out.print(i + " ");
			for (int j = 0; j < target[0].length; j++){
				for (int k = 0; k < gaussPts.length; k++)
					for (int m = 0; m < gaussPts[0].length; m++)
					{
						d = Distance.distance(i-gaussPts[k][m][0], j-gaussPts[k][m][1]);
						if (d < 5*spread){ 
							g = height*Math.exp(-(d*d)/(spread*spread));
							target[i][j] += g;
					}
				}
				if (limitTo1)
					target[i][j] = Math.min(target[i][j], 1);
			}
		}
	}
	public static void putGaussDots(double[][] target, double[][] gaussPts, double spread, double height)
	{
		double g;
		double d;
		FieldOps.zero(target);
		for (int i = 0; i < target.length; i++)
			for (int j = 0; j < target[0].length; j++){
				for (int k = 0; k < gaussPts.length; k++)
					{
						d = Distance.distance(i-gaussPts[k][0], j-gaussPts[k][1]);
						if (d < 5*spread){ 
							g = height*Math.exp(-(d*d)/(spread*spread));
							target[i][j] += g;
					}
				}
			}
	}
	
	public static void putGaussDots(double[][] target, AtomicCoordinatesSet latt, double spread, double height)
	{
		FieldOps.zero(target);
		int N = target.length;
		 double[] atomc1 = latt.getAtomicCoords(-N/2, -N/2);
		 double[] atomc2 = latt.getAtomicCoords(-N/2, 3*N/2);
		 double[] atomc3 = latt.getAtomicCoords(3*N/2, -N/2);
		 double[] atomc4 = latt.getAtomicCoords(3*N/2, 3*N/2);
		 double[] min = new double [2], max = new double [2];
		 min[0] = Math.min(Math.min(atomc1[0], atomc2[0]),Math.min(atomc3[0], atomc4[0]));
		 min[1] = Math.min(Math.min(atomc1[1], atomc2[1]),Math.min(atomc3[1], atomc4[1]));
		 max[0] = Math.max(Math.max(atomc1[0], atomc2[0]),Math.max(atomc3[0], atomc4[0]));
		 max[1] = Math.max(Math.max(atomc1[1], atomc2[1]),Math.max(atomc3[1], atomc4[1]));
		 
		 int xmn = FieldOps.round(min[0]);
		 int xmx = FieldOps.round(max[0])+1;
		 int ymn = FieldOps.round(min[1]);
		 int ymx = FieldOps.round(max[1])+1;
		 
		 for (int i = xmn; i < xmx; i++)
			 for (int j = ymn; j < ymx; j++)
			 {
				 double[] r = latt.getPixelCoords(i, j);
				 
				 for (int p = (int)(r[0]-5*spread); p < (int)(r[0]+5*spread); p++)
					 for (int q = (int)(r[1]-5*spread); q < (int)(r[1]+5*spread); q++)
						 if (FieldOps.withinBounds(p, q, target)){
							 double d = Distance.distance(p-r[0], q-r[1]);
							 target[p][q] += height*Math.exp(-(d*d)/(spread*spread));;
						 }
				 
			 }
		 
//			for (int j = 0; j < target[0].length; j++){
//				for (int k = 0; k < gaussPts.length; k++)
//					{
//						d = Distance.distance(i-gaussPts[k][0], j-gaussPts[k][1]);
//						if (d < 5*spread){ 
//							g = height*Math.exp(-(d*d)/(spread*spread));
//							target[i][j] += g;
//					}
//				}
			
	}
	
	//this is to generate the field which will shift a pixel with certain atomic coordinates in "imperfect" to the place cooresponding to the same 
	//atomic coordinate values in "perfect." For use in DriftCorrectionConsole
	public static void putUField(AtomicCoordinatesSet imperfect, AtomicCoordinatesSet perfect, double[][][] u)
	{
		double[] coordImp = new double [2];
		double[] pixPerf = new double [2];
//		String[] log = new String[u.length*u[0].length];
		
		for (int i = 0; i < u.length; i++)
			for (int j = 0; j < u[0].length; j++)
			{
				coordImp = imperfect.getAtomicCoords(i, j);
				pixPerf = perfect.getPixelCoords(coordImp);
//				if (pixPerf[0] == Double.NaN || pixPerf[1] == Double.NaN)
//				{
//					log[i*u.length + j] = "" + i +  "\t" + j + "\t" + imperfect.forDisplay2();
//				}
				u[i][j][0] = (pixPerf[0] - i);
				u[i][j][1] = (pixPerf[1] - j);
//				if (Complex.mag(u[i][j]) > 10)
//				{
//					System.out.println("wtf?");
//				}
			}
//		FileOps.writeLines(log);
	}
	
	public static void add(double[][][] a, double[][][] b, double[][][] target)
	{
		for (int i = 0; i < a.length; i++)
			for (int j = 0; j < a[i].length; j++)
				for (int k = 0; k < a[i][j].length; k++)
					target[i][j][k] = a[i][j][k] + b[i][j][k];
	}
	public static void addVectorTo(double[][][] a, double[] vector)
	{
		for (int i = 0; i < a.length; i++)
			for (int j = 0; j < a[i].length; j++)
				for (int k = 0; k < a[i][j].length; k++)
					a[i][j][k] += vector[k];
	}
	public static void add(double[][] a, double[][] b, double[][] target)
	{
		if (a == null)
			a = new double [b.length][b[0].length];
		for (int i = 0; i < a.length; i++)
			for (int j = 0; j < a[i].length; j++)
					target[i][j] = a[i][j] + b[i][j];
	}
	public static double[][] add(double[][] a, double[][] b)
	{
		double[][] target = new double [a.length][a[0].length];
		for (int i = 0; i < a.length; i++)
			for (int j = 0; j < a[i].length; j++)
					target[i][j] = a[i][j] + b[i][j];
		return target;
	}
	
	//sums up all the layers in data into target
	public static void sumUp(double[][][] data, double[][] target)
	{
		for (int j = 0; j < data[0].length; j++)
			for (int k = 0; k < data[0][j].length; k++)
				target[j][k] = 0;
		
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[i].length; j++)
				for (int k = 0; k < data[i][j].length; k++)
					target[j][k] += data[i][j][k];
	}

	public static double[][][] splitByLineHorizontal(double[][] data) {
		double[][][] ans = new double [2][data.length][data[0].length];
			for (int i = 0; i < data.length; i++)
				for (int j = 0; j < data[i].length; j++)
					ans[j%2][i][j] = data[i][j];
		FieldOps.changeZeroToAverage(ans[0]);
		FieldOps.changeZeroToAverage(ans[1]);
		return ans;
	}

	public static double[][] transpose(double[][] table) {
		double[][] ans = new double [table[0].length][table.length];
		for (int i = 0; i < ans.length; i++)
			for (int j= 0; j < ans[0].length; j++)
			{
				ans[i][j] = table[j][i];
			}
		return ans;
	}
	public static int[][] transpose(int[][] table) {
		int[][] ans = new int [table[0].length][table.length];
		for (int i = 0; i < ans.length; i++)
			for (int j= 0; j < ans[0].length; j++)
			{
				ans[i][j] = table[j][i];
			}
		return ans;
	}
	public static void putTranspose(double[][] table) {
		double[][] ans = new double [table[0].length][table.length];
		for (int i = 0; i < ans.length; i++)
			for (int j= 0; j < ans[0].length; j++)
			{
				ans[i][j] = table[j][i];
			}
		for (int i = 0; i < ans.length; i++)
			for (int j= 0; j < ans[0].length; j++)
			{
				table[i][j] = ans[i][j];
			}
	}
	public static double[][] rotateMinus90(double[][] table) {
		double[][] ans = new double [table[0].length][table.length];
		for (int i = 0; i < ans.length; i++)
			for (int j= 0; j < ans[0].length; j++)
			{
				ans[i][j] = table[j][ans.length-i-1];
			}
		return ans;
	}
	public static double[][] rotatePlus90(double[][] table) {
		double[][] ans = new double [table[0].length][table.length];
		for (int i = 0; i < ans.length; i++)
			for (int j= 0; j < ans[0].length; j++)
			{
				ans[i][j] = table[ans[0].length-j-1][i];
			}
		return ans;
	}
	/**
	 * rotates about +90 keeping the selected pixel fixed.
	 * @param table
	 * @return
	 */
	public static double[][] rotatePlus90_aboutPixel(double[][] table, int cx, int cy) {
		double[][] ans = new double [table[0].length][table.length];
		int x, y, xp, yp, ip, jp, nx = table.length, ny = table[0].length;
		for (int i = 0; i < ans.length; i++)
			for (int j= 0; j < ans[0].length; j++)
			{
				x = i-cx; y = j-cy;
				yp = -x;
				xp = y;
				ip = ((xp + cx) + nx)%nx;
				jp = ((yp + cy) + ny)%ny;
				ans[i][j] = table[ip][jp];
			}
		return ans;
	}
	
	/**
	 * This fits the z[i][j] to the form a + b*i + c*j and returns the array {a, b, c}.
	 * The fitting is based on the three linear equations
	 * D = N^2 a + N Sx b + N Sy c
	 * Dx = N Sx a + N Vx b + Vxy c
	 * Dy = N Sy a + Vxy b + N Vy c
	 * 
	 * which are based on minimizing R^2 = Sum(i, j) [z[i][j] - (a + b*i + c*j)]^2 with respect to a, b, and c
	 * (compare with http://mathworld.wolfram.com/LeastSquaresFitting.html).
	 * 
	 * Here D = Sum(i, j) z[i][j], Dx = Sumij(z*x[i]), Dy = Sumij(z*y[i])
	 * Sx = Sumi(x[i]), Sy = Sumj(y[j]), Vx = Sumi(x[i]^2), Vy = Sumj(y[j])^2, Vxy = Sumij(x[i]*y[j])
	 * N = the number (array is assumed square).
	 * 
	 * Throughout we assume that x[i] = i and y[j] = j
	 * 
	 * @param z
	 * @return
	 */
	public static double[] fitToPlaneSquare(double[][] z)
	{
		double sx = 0, sy =0, vx = 0, vy = 0, vxy = 0;
		double d = 0, dx = 0, dy = 0;
		int N = z.length;
		int M = z[0].length;
		
		double a, b, c;
		//If N != M next for loop must be split.
		for (int i = 0; i < N; i++)
		{
			sx += i;
			vx += i*i;
		}
		for (int j = 0; j < M; j++)
		{
			sy += j;
			vy += j*j;
		}
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
			{
				d += z[i][j];
				dx += z[i][j]*i;
				dy += z[i][j]*j;
				vxy += i*j;
			}
		
		double n1x = (sx/N - vx/sx);
		double n1y = (sy/M - vy/sy);
		double n2x = (d/(N*N) - dx/(N*sx));
		double n2y = (d/(N*N) - dy/(M*sy));
//		double n31 = (vxy/(N*sx) - sy/N); //both n31 and n32 are identically 0.
//		double n32 = (vxy/(N*sy) - sx/N);
		b = n2x/n1x;
		c = n2y/n1y;
		a = d/(N*N) - (sx/N)*b - (sy/M)*c;
		System.out.println(a + "\t" + b + "\t" + c);
		return new double[] {a, b, c};
	}
	/**
	 * This method treats the double-boolean combination as a cloud of points in which 
	 * only the points [i][j] for which the boolean is true are included. If the boolean
	 * is null, everything is included.
	 * Code adapted from
	 * http://stackoverflow.com/questions/9243645/weighted-least-square-fit-a-plane-to-3d-point-set
	 * @param z
	 * @return
	 */
	public static double[] fitToPlane(double[][] z, boolean[][] include)
	{
		if (include == null){
			include = new boolean[z.length][z[0].length];
			for (int i = 0; i < z.length; i++)
				for (int j = 0; j < z[0].length; j++)
					include[i][j] = true;
		}

		double sx = 0, sy = 0, sz = 0;
		double sxx = 0, sxy = 0, sxz = 0, syy = 0, syz = 0;
		int n = 0;
		for (int i = 0; i < z.length; i++)
			for (int j = 0; j < z[0].length; j++)
			{
				if (include[i][j])
				{
					sx += i;
					sy += j;
					sz += z[i][j];
					sxx += i*i;
					syy += j*j;
					sxy += j*i;
					sxz += i*z[i][j];
					syz += j*z[i][j];
					n++;
				}
			}
		double[][] A = {
				{sxx, sxy, sx},
				{sxy, syy, sy},
				{sx, sy, n}};
		double[] B = {sxz,syz,sz};		
		//use Apahce linear solver:
		RealMatrix ac = new Array2DRowRealMatrix(A, false);
		DecompositionSolver solver = new LUDecomposition(ac).getSolver();
		RealVector bs = new ArrayRealVector(B, false);

		RealVector solution = null;
		try{solution = solver.solve(bs);
		}
		catch (Exception e)
		{
			return new double [] {sz/n, 0, 0};
		}
		double[] answer = solution.toArray();
//		Printer.printlnHorizontal(answer);
		return new double[] {answer[2],answer[0],answer[1]};
	}
	public static double[] fitToPlane(double[] x, double[] y, double[] z)
	{

		double sx = 0, sy = 0, sz = 0;
		double sxx = 0, sxy = 0, sxz = 0, syy = 0, syz = 0;
		int n = 0;
		for (int i = 0; i < z.length; i++)
			{
					sx += x[i];
					sy += y[i];
					sz += z[i];
					sxx += x[i]*x[i];
					syy += y[i]*y[i];
					sxy += y[i]*x[i];
					sxz += x[i]*z[i];
					syz += y[i]*z[i];
					n++;
			}
		double[][] A = {
				{sxx, sxy, sx},
				{sxy, syy, sy},
				{sx, sy, n}};
		double[] B = {sxz,syz,sz};		
		//use Apahce linear solver:
		RealMatrix ac = new Array2DRowRealMatrix(A, false);
		DecompositionSolver solver = new LUDecomposition(ac).getSolver();
		RealVector bs = new ArrayRealVector(B, false);

		RealVector solution = null;
		try{solution = solver.solve(bs);
		}
		catch (Exception e)
		{
			return new double [] {sz/n, 0, 0};
		}
		double[] answer = solution.toArray();
//		Printer.printlnHorizontal(answer);
		return new double[] {answer[2],answer[0],answer[1]};
	}
	/**
	 * This code is based on math found at 
	 * http://www.efunda.com/math/leastsquares/lstsqr2dcurve.cfm
	 * adapted to the 2D case.
	 * 
	 * 
	 * @param z
	 * @param include
	 * @return
	 */
	public static double[] fitToParabola(double[][] z, boolean[][] include)
	{
		if (include == null){
			include = new boolean[z.length][z[0].length];
			for (int i = 0; i < z.length; i++)
				for (int j = 0; j < z[0].length; j++)
					include[i][j] = true;
		}

		double sx = 0, sy = 0, sz = 0;
		double sxx = 0, sxy = 0, sxz = 0, syy = 0, syz = 0;
		double sxxx = 0, sxxy = 0, sxyy= 0, syyy = 0;
		double sxxz = 0, sxyz = 0, syyz = 0;
		double sxxxx = 0, sxxxy = 0, sxxyy = 0, sxyyy = 0, syyyy = 0;
		int n = 0;
		double x, y;
		for (int i = 0; i < z.length; i++)
			for (int j = 0; j < z[0].length; j++)
			{
				if (include[i][j])
				{	
					x = (double)i;
					y = (double)j;
					sx += x;
					sy += y;
					sz += z[i][j];
					sxx += x*x;
					syy += y*y;
					sxy += y*x;
					sxz += x*z[i][j];
					syz += y*z[i][j];
					sxxx += x*x*x;
					sxxy += x*x*y;
					sxyy += x*y*y;
					syyy += y*y*y;
					sxxxx += x*x*x*x;
					sxxxy += x*x*x*y;
					sxxyy += x*x*y*y;
					sxyyy += x*y*y*y;
					syyyy += y*y*y*y;
					sxyz += x*y*z[i][j];
					sxxz += x*x*z[i][j];
					syyz += y*y*z[i][j];
					n++;
				}
			}
		double[][] A = {
				{n, sx, sy, sxy, sxx, syy},
				{sx, sxx, sxy, sxxy, sxxx, sxyy},
				{sy, sxy, syy, sxyy, sxxy, syyy},
				{sxy, sxxy, sxyy, sxxyy, sxxxy, sxyyy},
				{sxx, sxxx, sxxy, sxxxy, sxxxx, sxxyy},
				{syy, sxyy, syyy, sxyyy, sxxyy, syyyy}};
//		System.out.println(Printer.getTable(A));
//		A = FieldOps.transpose(A);
		double[] B = {sz, sxz, syz, sxyz, sxxz, syyz};		
		//use Apahce linear solver:
		RealMatrix ac = new Array2DRowRealMatrix(A, false);
		DecompositionSolver solver = new LUDecomposition(ac).getSolver();
		RealVector bs = new ArrayRealVector(B, false);

		RealVector solution = null;
		try{solution = solver.solve(bs);
		}
		catch (Exception e)
		{
			return fitToPlane(z, include);
		}
		double[] answer = solution.toArray();
//		Printer.printlnHorizontal(answer);
		return answer;
	}
	public static double[] fitToNSheet(double[][] z, boolean[][] include, int power)
	{
		if (include == null){
			include = new boolean[z.length][z[0].length];
			for (int i = 0; i < z.length; i++)
				for (int j = 0; j < z[0].length; j++)
					include[i][j] = true;
		}

		double[][] paramsNoZ = new double [2*power+1][];
		double[][] paramsZ = new double [power+1][];
		
		double[][] dparamsNoZ = new double [2*power+1][];
		double[][] dparamsZ = new double [power+1][];
		//The paramsNoZ table is x, y (i=0), (xx, xy, yy) (i = 1), etc.
		//The paramsZ table is z (i=0), xz, yz, (i=1), etc.
		for (int i = 0; i < power+1; i++)
		{
			paramsZ[i] = new double [i+1];
			dparamsZ[i] = new double [i+1];
		}
		for (int i = 0; i < 2*power+1; i++)
		{
			paramsNoZ[i] = new double [i+1];
			dparamsNoZ[i] = new double [i+1];
		}
		
		double x, y;
		for (int i = 0; i < z.length; i++)
			for (int j = 0; j < z[0].length; j++)
				if (include[i][j])
				{	
					x = (double)(i-z.length/2);
					y = (double)(j-z[0].length/2);
					for (int p = 0; p < 2*power+1; p++){
						for (int q = 0; q < dparamsNoZ[p].length; q++)
							dparamsNoZ[p][q] = Math.pow(x, dparamsNoZ[p].length-(q+1)) * Math.pow(y, q);

						for (int q = 0; q < dparamsNoZ[p].length; q++)
							paramsNoZ[p][q] += dparamsNoZ[p][q];
					}
					for (int p = 0; p < power+1; p++){
						for (int q = 0; q < dparamsZ[p].length; q++)
							dparamsZ[p][q] = Math.pow(x, dparamsZ[p].length-(q+1)) * Math.pow(y, q) * z[i][j];
						for (int q = 0; q < dparamsZ[p].length; q++)
							paramsZ[p][q] += dparamsZ[p][q];
					}
				}
				
		double[][] matrix = new double [((power+1)*(power+2))/2][((power+1)*(power+2))/2];
		
		
//		double[] paramsNoZTable = FieldOps.getArray(paramsNoZ);
		int p0 = 0, q0 = 0;
		int dp = 0, dq = 0;
		for (int i = 0; i < matrix.length; i++){
			p0 = (int)((-1 + Math.sqrt(1 + 8*i))/2);
			q0 = i - (p0*(p0+1))/2;
			
			for (int j = 0; j < matrix[0].length; j++)
			{
				dp = (int)((-1 + Math.sqrt(1 + 8*j))/2);
				dq = j - (dp*(dp+1))/2;
				matrix[i][j] = paramsNoZ[p0+dp][q0+dq]; 
			}
			
		}
//		System.out.println(Printer.getTable(A));
//		A = FieldOps.transpose(A);
		double[] B = new double [((power+1)*(power+2))/2];
		for (int i = 0; i < matrix.length; i++)
		{
			p0 = (int)((-1 + Math.sqrt(1 + 8*i))/2);
			q0 = i - (p0*(p0+1))/2;
			B[i] = paramsZ[p0][q0];
		}
		//use Apahce linear solver:
		RealMatrix ac = new Array2DRowRealMatrix(matrix, false);
		DecompositionSolver solver = new LUDecomposition(ac).getSolver();
		RealVector bs = new ArrayRealVector(B, false);

		RealVector solution = null;
		try{solution = solver.solve(bs);
		}
		catch (Exception e)
		{
			return fitToPlane(z, include);
		}
		double[] answer = solution.toArray();
//		Printer.printlnHorizontal(answer);
		return answer;
	}
	/**
	 * This assumes that the layer is in chronological order i.e. no rotation of the lines or anything.
	 */
	public static void subtractLinePolynomialFit(double[][] z, int power, boolean insertSpaces)
	{
		int nx = z.length, ny = z[0].length;
		double[] y = new double [nx*ny];
		double[] x = new double [nx*ny];
		int n = 0;
		for (int j = 0; j < ny; j++)
			for (int i = 0; i < nx; i++)
			{
				y[n] = z[i][j];
				x[n] = insertSpaces ? n + j*nx : n;
				n++;
			}
		
		double[] ans = ArrayOps.subtractPolynomialFit(x, y, power)[0];
		n = 0;
		for (int j = 0; j < ny; j++)
			for (int i = 0; i < nx; i++)
				z[i][j] = ans[n++];
	}
	public static void subtractPlaneFit(double[][] z)
	{
		double[] a = fitToPlane(z, null);
		for (int i = 0; i < z.length; i++)
			for (int j = 0; j < z[0].length; j++)
			{
				z[i][j] -= 0.9999*a[0] + i*a[1] + j*a[2];
			}
	}
	public static void subtractParabolicFit(double[][] z)
	{
		double[] a = fitToParabola(z, null);
		for (int i = 0; i < z.length; i++)
			for (int j = 0; j < z[0].length; j++)
			{
				z[i][j] -= 0.9999*a[0] + i*a[1] + j*a[2] + i*j*a[3] + i*i*a[4] + j*j*a[5];
			}
	}
	public static void subtractNSheetFit(double[][] z, int power)
	{
		double[] a = fitToNSheet(z, null, power);
		a[0] = 0.9999*a[0];
		int n = 0;
		int nxhalf = z.length/2, nyhalf = z[0].length/2;
		for (int i = 0; i < z.length; i++)
			for (int j = 0; j < z[0].length; j++)
			{
				n = 0;
				for (int p = 0; p < power+1; p++){
					for (int q = 0; q < p+1; q++)
						z[i][j] -= Math.pow((double)(i-nxhalf), (p - q)) * Math.pow((double)(j-nyhalf), q) * a[n++];
				}
			}
		}

	public static void subtractNSheetFit(double[][] z, double[] a, int power)
	{
		a[0] = 0.9999*a[0];
		int n = 0;
		int nxhalf = z.length/2, nyhalf = z[0].length/2;
		for (int i = 0; i < z.length; i++)
			for (int j = 0; j < z[0].length; j++)
			{
				n = 0;
				for (int p = 0; p < power+1; p++){
					for (int q = 0; q < p+1; q++)
						z[i][j] -= Math.pow((double)(i-nxhalf), (p - q)) * Math.pow((double)(j-nyhalf), q) * a[n++];
				}
			}
		}
	public static void putNSheetFit(double[][] z, double[] a, int power)
	{
		a[0] = 0.9999*a[0];
		int n = 0;
		int nxhalf = z.length/2, nyhalf = z[0].length/2;
		for (int i = 0; i < z.length; i++)
			for (int j = 0; j < z[0].length; j++)
			{
				z[i][j] = 0;
				n = 0;
				for (int p = 0; p < power+1; p++){
					for (int q = 0; q < p+1; q++)
						z[i][j] += Math.pow((double)(i-nxhalf), (p - q)) * Math.pow((double)(j-nyhalf), q) * a[n++];
				}
			}
		}
	/**
	 * This subtracts a different plane fit from each block of pixels
	 * @param z
	 * @param blocks
	 */
	public static void subtractPlaneFits(double[][] z, int[][] blocks)
	{
		int max = FieldOps.max(blocks);
		int min = FieldOps.min(blocks);
		
		boolean[][] include = new boolean[z.length][z[0].length];
		for (int k = min; k <= max; k++)
		{
			for (int i = 0; i < z.length; i++)
				for (int j = 0; j < z[0].length; j++)
					include[i][j] = blocks[i][j] == k;
			
		double[] a = fitToPlane(z, include);
		for (int i = 0; i < z.length; i++)
			for (int j = 0; j < z[0].length; j++)
				if (include[i][j])
					z[i][j] -= 0.9999*a[0] + i*a[1] + j*a[2];
		}
		
	}
	
	/**
	 * Assumes data and mask are the SAME SIZE.
	 * @param data
	 * @param mask
	 * @return
	 */
	public static double[][] convolveFourier(double[][] data, double[][] mask)
	{
		int nx = data.length, ny = data[0].length; 
		double[][][] fftzData = FFTOps.obtainFFT(data).fHat;
		double[][][] fftzMask = FFTOps.obtainFFT(mask).fHat;
		double[][][] ansFFTZ = new double [nx][ny][2];
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				Complex.product(fftzData[i][j], fftzMask[i][j], ansFFTZ[i][j]);
			}
		double[][][] ansTemp = new double [nx][ny][2];
		FFTOps.putIFFT(ansFFTZ, ansTemp, false);
		return getIndex(ansTemp, 0);
	}
	public static double[][] getConvolution(double[][] data, double[][] mask)
	{
		return convolveFourier(data, mask);
	}
	
	public static void subtractPlaneFitsNoMean(double[][] z, int[][] blocks)
	{
		int max = FieldOps.max(blocks);
		int min = FieldOps.min(blocks);
		
		boolean[][] include = new boolean[z.length][z[0].length];
		for (int k = min; k <= max; k++)
		{
			for (int i = 0; i < z.length; i++)
				for (int j = 0; j < z[0].length; j++)
					include[i][j] = blocks[i][j] == k;
			
		double[] a = fitToPlane(z, include);
		for (int i = 0; i < z.length; i++)
			for (int j = 0; j < z[0].length; j++)
				if (include[i][j])
					z[i][j] -= i*a[1] + j*a[2];
		}
		
	}
	public static void subtractParabolaFits(double[][] z, int[][] blocks)
	{
		int max = FieldOps.max(blocks);
		int min = FieldOps.min(blocks);
		
		boolean[][] include = new boolean[z.length][z[0].length];
		for (int k = min; k <= max; k++)
		{
			for (int i = 0; i < z.length; i++)
				for (int j = 0; j < z[0].length; j++)
					include[i][j] = blocks[i][j] == k;
			
		double[] a = fitToParabola(z, include);
		for (int i = 0; i < z.length; i++)
			for (int j = 0; j < z[0].length; j++)
				if (include[i][j])
					z[i][j] -= 0.9999*a[0] + i*a[1] + j*a[2] + a[3]*i*j + a[4]*i*i + a[5]*j*j;
		}
		
	}
	public static double getSpectralFunction(double w, double re, double im, double bareEnergy)
	{
		return (im)/(Math.pow(w - bareEnergy - re, 2) + im*im);
	}

	public static void multiply(double[][][] z, double f) {
		// TODO Auto-generated method stub
		for (int i = 0; i < z.length; i++)
			for (int j = 0; j < z[0].length; j++)
			{
				z[i][j][0] *= f;
				z[i][j][1] *= f;
			}
	}
	
	/**
	 * This can be used for suppressing FFT modes.
	 * @param z
	 * @param f
	 */
	public static void multiply(double[][][] z, double[][] f)
	{
		for (int i = 0; i < z.length; i++)
			for (int j = 0; j < z[0].length; j++)
			{
				z[i][j][0] *= f[i][j];
				z[i][j][1] *= f[i][j];
			}
		
	}
	public static void multiply(double[][][] z, double[][] f, double[][][] target)
	{
		for (int i = 0; i < z.length; i++)
			for (int j = 0; j < z[0].length; j++)
			{
				target[i][j][0] = z[i][j][0] * f[i][j];
				target[i][j][1] = z[i][j][1] * f[i][j];
			}
		
	}
	public static boolean withinBounds(int i, int j, double[][] array)
	{
		return i >= 0 && j >= 0 && i < array.length && j < array[i].length;
	}
	public static boolean withinBounds(int i, int j, boolean[][] array)
	{
		return i >= 0 && j >= 0 && i < array.length && j < array[i].length;
	}
	public static void replaceWithGradMag(double[][] data) {
		// TODO Auto-generated method stub
		double[][][] grad = FieldOps.gradientNM2(data);
		double[][] mag = FieldOps.magnitude(grad);
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
				data[i][j] = mag[i][j];
	}
	public static double[][] gradMag(double[][] data) {
		// TODO Auto-generated method stub
		double[][][] grad = FieldOps.gradientNM2(data);
		return FieldOps.magnitude(grad);

		
	}
	public static double[] getArray(double[][] data)
	{
		int n = 0;
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[i].length; j++)
				n++;
		double[] ans = new double [n];
		
		n = 0;
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[i].length; j++)
				ans[n++] = data[i][j];
		return ans;
	}
	public static double[] getArray(double[][][] data)
	{
		int n = 0;
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[i].length; j++)
				for (int k = 0; k < data[i][j].length; k++)
					n++;
		double[] ans = new double [n];
		
		n = 0;
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[i].length; j++)
				for (int k = 0; k < data[i][j].length; k++)
					ans[n++] = data[i][j][k];
		return ans;
	}
	public static double[] getArray(double[][][] data, int imin, int imax)
	{
		int n = 0;
		for (int i = imin; i < imax; i++)
			for (int j = 0; j < data[i].length; j++)
				for (int k = 0; k < data[i][j].length; k++)
					n++;
		double[] ans = new double [n];
		
		n = 0;
		for (int i = imin; i < imax; i++)
			for (int j = 0; j < data[i].length; j++)
				for (int k = 0; k < data[i][j].length; k++)
					ans[n++] = data[i][j][k];
		return ans;
	}
	public static double[] getArray(double[][] data, boolean[][] usePixel)
	{
		int n = 0;
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[i].length; j++)
				if (usePixel[i][j])
					n++;
		double[] ans = new double [n];
		
		n = 0;
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[i].length; j++)
				if (usePixel[i][j])
					ans[n++] = data[i][j];
		return ans;
	}
	
	/**
	 * This replaces the lines by the FFT of the lines
	 * @return
	 */
	public static double[][] getReplacedRowFFT(double[][] data, boolean horizontal, boolean log)
	{
		DoubleFFT_1D fft = null;
		double[] temp = null;
		double[][] ans = new double [data.length][data[0].length];
		if (horizontal)
		{
			fft = new DoubleFFT_1D(data.length);
			temp = new double [2*data.length];
			for (int j = 0; j < data[0].length; j++)
			{
				for (int i = 0; i < data.length; i++)
				{
					temp[2*i] = data[i][j];
					temp[2*i+1] = 0;
				}
				fft.complexForward(temp);
				for (int i = 0; i < data.length; i++)
				{
					ans[i][j] = Complex.mag(temp[2*i], temp[2*i+1]);
				}
			}
			if (log)
				FieldOps.log(ans);
			return ans;
		}
		else
		{
			fft = new DoubleFFT_1D(data[0].length);
			temp = new double [2*data[0].length];
			for (int i = 0; i < data.length; i++)
			{
				for (int j = 0; j < data[0].length; j++)
				{
					temp[2*j] = data[i][j];
					temp[2*j+1] = 0;
				}
				fft.complexForward(temp);
				for (int j = 0; j < data.length; j++)
				{
					ans[i][j] = Complex.mag(temp[2*j], temp[2*j+1]);
				}
			}
			if (log)
				FieldOps.log(ans);
			return ans;
		}
	}
	
	/**
	 * This returns an int array indicating the position of all pixels in data with respect to the ascending-order
	 * array bounds.
	 * If data[i][j] < bounds[0], 0, if data[i][j] > bounds[0] but < bounds[1] 1, etc.
	 * 
	 * @param data
	 * @param bounds
	 * @return
	 */
	public static int[][] getPositionInArray(double[][] data, double[] bounds)
	{
		int[][] ans = new int[data.length][data[0].length];
		for (int i = 0; i < ans.length; i++)
			for (int j = 0; j < ans[0].length; j++)
			{
				ans[i][j] = -1;
				if (data[i][j] < bounds[0]) ans[i][j] = 0;
				for (int k = 1; k < bounds.length; k++)
					if (data[i][j] < bounds[k] && data[i][j] >= bounds[k-1]) ans[i][j] = k;
				if (ans[i][j] == -1)
					ans[i][j] = bounds.length;
			}
		return ans;
	}
	/**
	 * returns a boolean map indicating whether each pixel belongs to a contiguous FALSE
	 * blob in the visited field, the blob starting at (i, j).
	 * It is assumed that visited[i][j] is itself false.
	 * The visited array IS modified in this method.
	 * @param visited
	 * @return
	 */
	public static boolean[][] getContiguousBlob(int i, int j, boolean[][] visited)
	{
		ArrayList<int[]> blob = getContiguousBlobList(i, j, visited);
		
		boolean[][] ans = new boolean [visited.length][visited[0].length];
		for (int k = 0; k < blob.size(); k++)
			ans[blob.get(k)[0]][blob.get(k)[1]] = true;
		return ans;
	}
	public static boolean[][] getTrueArray(ArrayList<int[]> blob, int nx, int ny)
	{	
		boolean[][] ans = new boolean [nx][ny];
		for (int k = 0; k < blob.size(); k++)
			ans[blob.get(k)[0]][blob.get(k)[1]] = true;
		return ans;
	}
	public static boolean[][] getTrueArrayBlobs(ArrayList<ArrayList<int[]>> blobs, int nx, int ny)
	{	
		boolean[][] ans = new boolean [nx][ny];
		for (int k = 0; k < blobs.size(); k++)
			for (int j = 0; j < blobs.get(k).size(); j++)
			ans[blobs.get(k).get(j)[0]][blobs.get(k).get(j)[1]] = true;
		return ans;
	}
	public static int[][] getContiguousBlobsInt(ArrayList<ArrayList<int[]>> blobs, int nx, int ny)
	{	
		int[][] ans = new int [nx][ny];
		for (int k = 0; k < blobs.size(); k++)
			for (int j = 0; j < blobs.get(k).size(); j++)
				ans[blobs.get(k).get(j)[0]][blobs.get(k).get(j)[1]] = k+1;
		return ans;
	}
	public static ArrayList<int[]> getContiguousBlobList(int i, int j, boolean[][] visited)
	{
		ArrayList<int[]> queue = new ArrayList<int[]>();
		ArrayList<int[]> blob = new ArrayList<int[]>();
		int count = 0;
		int nx = visited.length, ny = visited[0].length;
		
		int[] currentPoint;
		int x, y;
		count++;
		blob = new ArrayList<int[]>();
		if (!visited[i][j]){
			queue.add(new int[] {i,j});
			blob.add(new int[] {i,j});
		}
		while(queue.size() > 0)
		{
			currentPoint = queue.remove(0);
			x = currentPoint[0]; y = currentPoint[1];
			visited[x][y] = true;
			if (x > 0 && !visited[x-1][y]){
				queue.add(new int[] {x-1,y});
				blob.add(new int[] {x-1,y});
				visited[x-1][y] = true;
			}
			if (y > 0 && !visited[x][y-1]){
				queue.add(new int[] {x,y-1});
				blob.add(new int[] {x,y-1});
				visited[x][y-1] = true;
			}
			if (x < nx-1 && !visited[x+1][y]){
				queue.add(new int[] {x+1,y});
				blob.add(new int[] {x+1,y});
				visited[x+1][y] = true;
			}
			if (y < ny-1 && !visited[x][y+1]){
				queue.add(new int[] {x,y+1});
				blob.add(new int[] {x,y+1});
				visited[x][y+1] = true;
			}
		}
		return blob;
	}
	
	/**
	 * This returns an array list containing all the pixels in a contiguous blob, which are of the 
	 * same sign as data[i][j];
	 * @param i
	 * @param j
	 * @param data
	 * @return
	 */
	public static ArrayList<int[]> getContiguousBlobList(int i, int j, double[][] data)
	{
		ArrayList<int[]> queue = new ArrayList<int[]>();
		ArrayList<int[]> blob = new ArrayList<int[]>();
		int count = 0;
		int nx = data.length, ny = data[0].length;
		
		int[] currentPoint;
		int x, y;
		count++;
		blob = new ArrayList<int[]>();
		queue.add(new int[] {i,j});
		blob.add(new int[] {i,j});
		boolean[][] visited = new boolean[nx][ny];
		int sign = data[i][j] > 0 ? 1 : -1;
		while(queue.size() > 0)
		{
//			System.out.println(queue.size());
			currentPoint = queue.remove(0);
			x = currentPoint[0]; y = currentPoint[1];
			visited[x][y] = true;
			if (x > 0 && !visited[x-1][y] && data[x-1][y] * sign > 0){
				queue.add(new int[] {x-1,y});
				blob.add(new int[] {x-1,y});
				visited[x-1][y] = true;
			}
			if (y > 0 && !visited[x][y-1] && data[x][y-1] * sign > 0){
				queue.add(new int[] {x,y-1});
				blob.add(new int[] {x,y-1});
				visited[x][y-1] = true;
			}
			if (x < nx-1 && !visited[x+1][y] && data[x+1][y] * sign > 0){
				queue.add(new int[] {x+1,y});
				blob.add(new int[] {x+1,y});
				visited[x+1][y] = true;
			}
			if (y < ny-1 && !visited[x][y+1] && data[x][y+1] * sign > 0){
				queue.add(new int[] {x,y+1});
				blob.add(new int[] {x,y+1});
				visited[x][y+1] = true;
			}
		}
		return blob;
	}
	public static ArrayList<ArrayList<int[]>> getAllContiguousBlobs(boolean[][] visited)
	{
		ArrayList<ArrayList<int[]>> blobs = new ArrayList<ArrayList<int[]>>();
		ArrayList<int[]> blob;
		for (int i = 0; i < visited.length; i++)
			for (int j = 0; j < visited.length; j++){
				blob = getContiguousBlobList(i, j, visited);
				if (blob.size() > 0)
					blobs.add(blob);
			}
		return blobs;
	}
	public static ArrayList<int[]> getTruePoints(boolean[][] b)
	{
		ArrayList<int[]> ans = new ArrayList<int[]>();
		for (int i = 0; i < b.length; i++)
			for (int j = 0; j < b[0].length; j++)
				if (b[i][j]) ans.add(new int[] {i, j});
		return ans;
	}
	public static int getNumTrue(boolean[][] b)
	{
		int n = 0;
		for (int i = 0; i < b.length; i++)
			for (int j = 0; j < b[0].length; j++)
				if (b[i][j]) n++;
		return n;
	}
	
	/**
	 * This produces an int[][] array of blocky pixels. 
	 * @param isStepEdge
	 * @return
	 */
	public static int[][] splitByStepEdges(boolean[][] isStepEdge)
	{
		boolean[][] part;
		boolean[][] visited = isStepEdge.clone();
		int[][] ans = new int [visited.length][visited[0].length];
		for (int i = 0; i < visited.length; i++)
			for (int j = 0; j < visited[0].length; j++)
				ans[i][j] = -1;
		
		int k = 0;
		for (int i = 0; i < visited.length; i++)
			for (int j = 0; j < visited[0].length; j++)
				if (!visited[i][j])
				{
					part = getContiguousBlob(i, j, visited);
					for (int m = 0; m < visited.length; m++)
						for (int n = 0; n < visited[0].length; n++)
							if (part[m][n]) ans[m][n] = k; else;
					k++;
				}
		return ans;
	}
	
	/**
	 * Splits the field data into n fields, where n is the number of bins. Each data[i] contains all the data in the matrix such that bins[p][q] = min(bins)+i;
	 * @param data
	 * @param bins
	 * @return
	 */
	public static double[][][] splitByBins(double[][] data, int[][] bins)
	{
		int min = FieldOps.min(bins);
		int max = FieldOps.max(bins);
		int n = (max-min)+1;
		
		double[][][] ans = new double [n][data.length][data[0].length];
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
			{
				ans[bins[i][j]-min][i][j] = data[i][j];
			}
		return ans;
	}
	/**
	 * This returns an array of weight functions for the various bins. The weight functions are smoothed by the length scale L.
	 * @param data
	 * @param bins
	 * @return
	 */
	public static double[][][] getSmoothedWeightingFunctions(int[][] bins, double L)
	{
		int min = FieldOps.min(bins);
		int max = FieldOps.max(bins);
		int n = (max-min)+1;
		
		double[][][] ans = new double [n][bins.length][bins[0].length];
		for (int k = 0; k < n; k++){
			for (int i = 0; i < bins.length; i++)
				for (int j = 0; j < bins[0].length; j++)
				{
					ans[k][i][j] = (bins[i][j] == k) ? 1 : 0;
				}
			ans[k] = FieldOps.gaussSmooth(ans[k], L);
		}
		return ans;
	}
	public static double[][] getSmoothedWeightingFunctions(int[][] bins, double L, int index)
	{
		int min = FieldOps.min(bins);
		int max = FieldOps.max(bins);
		int n = (max-min)+1;
		
		double[][] ans = new double [bins.length][bins[0].length];
			for (int i = 0; i < bins.length; i++)
				for (int j = 0; j < bins[0].length; j++)
				{
					ans[i][j] = (bins[i][j] == index) ? 1 : 0;
				}
		ans = FieldOps.gaussSmooth(ans, L);
		return ans;
	}
	
	/**
	 * returns an array where ans[i][j] is the index of the bin field[i][j] falls into. indices go from 0 to nbins-1
	 * @param field
	 * @param nbins
	 * @return
	 */
	public static int[][] getPercentileBinsForField(double[][] field, int nbins)
	{
		int nx = field.length, ny = field[0].length;
		int[][] ans = new int [nx][ny];
		
		double[] array = FieldOps.getArray(field);
		ArrayOps.quicksort(array);
		if (nbins <= 1) return ans;
		
		double[] borders = new double [nbins-1];
		for (int i = 0; i < nbins-1; i++){
			borders[i] = array[(i+1)*array.length/(nbins)];
			System.out.println(borders[i]);
		}
		
		for (int i = 0 ; i < nx; i++)
			for (int j = 0; j < ny; j++)
				for (int k = 1; k < nbins; k++)
				{
					if (field[i][j] >= borders[k-1])
						ans[i][j] = k;
				}
		return ans;
	}
	/**
	 * This sorts the pixels where usePixels is true and splits that mass into nbins bins of equal size.
	 * In the return array the true pixels are labled from 0 to nbins-1 and the false pixels are labeled with -1.
	 * @param field
	 * @param nbins
	 * @param usePixel
	 * @return
	 */
	public static int[][] getPercentileBinsForField(double[][] field, int nbins, boolean[][] usePixel)
	{
		int nx = field.length, ny = field[0].length;
		int[][] ans = new int [nx][ny];
		
		double[] array = FieldOps.getArray(field, usePixel);
		ArrayOps.quicksort(array);
		if (nbins <= 1) return ans;
		
		double[] borders = new double [nbins-1];
		for (int i = 0; i < nbins-1; i++)
			borders[i] = array[(i+1)*array.length/(nbins)];
		
		for (int i = 0 ; i < nx; i++)
			for (int j = 0; j < ny; j++)
				if (usePixel[i][j])
					for (int k = 1; k < nbins; k++)
					{
						if (field[i][j] >= borders[k-1])
							ans[i][j] = k;
					}
				else
					ans[i][j] = -1;
		return ans;
	}
	public static double[] getPercentiles(double[][] field, double minperc, double maxperc)
	{
		double[] sorted = getSortedDump(field);
		int n = sorted.length-1;
		double min = sorted[(int)(Math.max(0, minperc)*n)];
		double max = sorted[(int)(Math.min(1, maxperc)*n)];
		
		return new double [] {min, max};
	}
	/**
	 * Resizes by the factor. It is assumed that the factor is an integer or the reciprocal of an integer.
	 * @param truePix
	 * @param imageOverField
	 * @return
	 */
	public static boolean[][] resizeByFactor(boolean[][] truePix, double imageOverField) {
		int scale = 1;
		if (imageOverField > 1)
			scale = (int)imageOverField;
		else if (imageOverField < 1)
			scale = FieldOps.round(1/imageOverField);
		
		if (scale == 1)
			return truePix.clone();
		else if (imageOverField < 1) //shrink the output down
		{
			boolean[][] ans = new boolean [truePix.length/scale][truePix[0].length/scale];
			for (int i = 0; i < ans.length; i++)
				for (int j = 0; j < ans[0].length; j++)
					ans[i][j] = truePix[scale*i][scale*j];
			return ans;
		}
		else
		{
			boolean[][] ans = new boolean [truePix.length*scale][truePix[0].length*scale];
			for (int i = 0; i < ans.length; i++)
				for (int j = 0; j < ans[0].length; j++)
					ans[i][j] = truePix[i/scale][j/scale];
			return ans;
			
		}
	}
	/**
	 * This adds a constant to all values in the dataset such that the center pixel of the dataset is equal to an integer multiple of mod.
	 * @param phaseCont
	 * @param mod
	 */
	public static void shiftCenterModulo(double[][] phaseCont, double mod) {
		double remainder = ((phaseCont[phaseCont.length/2][phaseCont[0].length/2] - mod/2) % mod) + mod/2;
		for (int i = 0; i < phaseCont.length; i++)
			for (int j = 0; j < phaseCont.length; j++)
				phaseCont[i][j] -= remainder;
	}
	
	public static BicubicSplineInterpolatingFunction getBicubicInterpoation(double[][] data)
	{
		double[] x = ArrayOps.generateArrayInclBoth(0, data.length-1, data.length);
		double[] y = ArrayOps.generateArrayInclBoth(0, data[0].length-1, data[0].length);
		double[][][] grad = FieldOps.gradient(data);
		double[][] cross = FieldOps.gradient(grad[0])[1];
		return new BicubicSplineInterpolatingFunction(x, y, data, grad[0], grad[1], cross);
	}
	public static double[] getCentroid(ArrayList<int[]> arrayList) {
		double x = 0, y = 0;
		for (int i = 0; i < arrayList.size(); i++)
		{
			x += arrayList.get(i)[0];
			y += arrayList.get(i)[1];
		}
		x /= arrayList.size();
		y /= arrayList.size();
		return new double[] {x, y};
	}
	/**
	 * We assume that in some sense f * g = h. We know h, and we know the function g, and we try to determine f.
	 * We use that hHat = fHat times gHat. Then fHat = hHat/gHat, and f is the inverse Fourier transform of fHat.
	 * Obviously this requires that gHat not be zero anywhere, and we will try to satisfy this. 
	 * @param h
	 * @param g
	 * @return
	 */
	public static double[][] getSimpleFourierDeconvolution(double[][] h, double[][] g){
		int nx = h.length, ny = h[0].length;
		double[][][] hHatZ = new double [nx][ny][2];
		double[][][] gHatZ = new double [nx][ny][2];
		FFTOps.putFFT(h, hHatZ, false);
		FFTOps.putFFT(g, gHatZ, false);
		double[][][] fHatZ = new double [nx][ny][2];
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				fHatZ[i][j] = Complex.getQuotient(hHatZ[i][j], gHatZ[i][j], 1000000);
			}
		double[][][] ans = new double [nx][ny][2];
		FFTOps.putIFFT(fHatZ, ans, false);
		return FieldOps.getIndex(ans, 0);
	}
	
	/**
	 * The expression for Wiener deconvolution is given in Wikipedia:
	 * G(f) = H*(f)S(f)/(H^2(f)S(f) + N(f))
	 * where the desired signal xHat is the convolution of g(t) with y(t).
	 * 
	 * The noise is assumed white so that N(f) = noiseLevel.
	 * The "mean power density of the input signal" S(f) is assumed
	 * to go to zero with the decay length GaussianDecay which removes resolution...
	 * @param h
	 * @param g
	 * @param noiseLevel
	 * @param gaussDecay
	 * @return
	 */
	public static double[][] getWienerDeconvolution(double[][] y, double[][] h, double noiseLevel, double gaussDecay, double noiseFloor){
		int nx = h.length, ny = h[0].length;
		double[][][] Hf = new double [nx][ny][2];
		double[][][] Yf = new double [nx][ny][2];
		FFTOps.putFFT(h, Hf, true);
		FFTOps.putFFT(y, Yf, true);
		double[][] HfMag2 = new double [nx][ny];
		double[][] Sf = new double [nx][ny];
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++){
				HfMag2[i][j] = Math.pow(Complex.mag(Hf[i][j]), 2);
				HfMag2[i][j] = Math.max(HfMag2[i][j], noiseFloor);
				double d = Distance.distance(i-nx/2, j-ny/2);
				Sf[i][j] = Math.pow(Complex.mag(Yf[i][j]), 2)*Math.exp(-d*d/(gaussDecay*gaussDecay));
			}
		double[][][] Gf = new double [nx][ny][2];
		double[][][] XHatf = new double [nx][ny][2];
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++){
				Gf[i][j][0] = (Hf[i][j][0])*Sf[i][j]/(HfMag2[i][j]*Sf[i][j] + noiseLevel);
				Gf[i][j][1] = (-Hf[i][j][1])*Sf[i][j]/(HfMag2[i][j]*Sf[i][j] + noiseLevel);
				XHatf[i][j] = Complex.product(Gf[i][j], Yf[i][j]);
			}
		double[][][] temp = new double [nx][ny][2];
		FFTOps.putIFFT(XHatf, temp, true);
		double[][] ans = FieldOps.getIndex(temp, 0);
		double[][] ansC = FieldOps.getIndex(temp, 0);
		FieldOps.shift(ansC, ans);
		return ans;
	}
	public static double[][][] getNM2Array(double[][][] ds) {
		double[][][] ans = new double [ds[0].length][ds[0][0].length][2];
		for (int i = 0; i < ans.length; i++)
			for (int j = 0; j < ans[0].length; j++)
			{
				ans[i][j][0] = ds[0][i][j];
				ans[i][j][1] = ds[1][i][j];
			}
		return ans;
	}
	public static void putNM2Array(double[][][] ds, double[][][] ans) {
		for (int i = 0; i < ans.length; i++)
			for (int j = 0; j < ans[0].length; j++)
			{
				ans[i][j][0] = ds[0][i][j];
				ans[i][j][1] = ds[1][i][j];
			}
	}
	public static double[][] spatialFilter(double[][] data, double[][] weight)
	{
		double mean = 0;
		int nx = data.length, ny = data[0].length;
		double[][] ans = new double[nx][ny];
		mean = FieldOps.mean(data, weight);
		System.out.println(mean);
		for (int j = 0; j < nx; j++)
			for (int k = 0; k < ny; k++)
			{
				ans[j][k] = weight[j][k]*data[j][k] + (1-weight[j][k])*mean;
			}
		return  ans;
	}
	
	public static double[][] spatialFilterZeros(double[][] data, double[][] weight)
	{
		double mean = 0;
		int nx = data.length, ny = data[0].length;
		double[][] ans = new double[nx][ny];
		System.out.println(mean);
		for (int j = 0; j < nx; j++)
			for (int k = 0; k < ny; k++)
			{
				ans[j][k] = weight[j][k]*data[j][k] + (1-weight[j][k])*mean;
			}
		return  ans;
	}
	
	public static void putAllTrues(boolean[][][] fftfilt, boolean[][] totalfilt) {
		for (int i = 0; i < totalfilt.length; i++)
			for (int j = 0; j < totalfilt.length; j++)
			{
				totalfilt[i][j] = false;
				for (int k = 0; k < fftfilt.length && !totalfilt[i][j]; k++)
					totalfilt[i][j] = totalfilt[i][j] || fftfilt[k][i][j];
			}
	}
	/**
	 * This returns the length scale where exp(-distance^2/L^2) has the average value desired percent.
	 * @param distance
	 * @param initL
	 * @param niter
	 * @param desiredPercent
	 * @return
	 */
	public static double getCorrectGaussianLength(double[][] distance, int niter, double desiredPercent){
		int n = 0;
		double lowerL = 0;
		double upperL = 4*distance.length;
		double L = (upperL+lowerL)/2;
		int nx = distance.length, ny = distance[0].length;
		double[][] gauss = new double [nx][ny];
		while (n <= niter){
			n++;
			for (int i = 0 ; i < nx; i++)
				for (int j = 0; j < ny; j++)
					gauss[i][j] = Math.exp(-(distance[i][j]*distance[i][j])/(L*L));
			double mean = FieldOps.mean(gauss);
			
			if (mean > desiredPercent){
				upperL = L;
			}
			else{
				lowerL = L;
			}
			L = (upperL+lowerL)/2;

			System.out.println(L);
		}
		return L;
	}
	public static double[][] getQuadratureSum(double[][] x, double[][] y) {
		double[][] d = new double [x.length][x[0].length];
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[i].length; j++)
				d[i][j] = Distance.distance(x[i][j], y[i][j]);
		return d;
	}

	public static double[][] getGaussianWeights(double[][] distance, double L){
		double[][] gauss = new double [distance.length][distance[0].length];
		for (int i = 0 ; i < distance.length; i++)
			for (int j = 0; j < distance[i].length; j++)
				gauss[i][j] = Math.exp(-(distance[i][j]*distance[i][j])/(L*L));
		return gauss;
	}
	public static class FieldWithDouble implements Comparable{

		public double v;
		public double[][] f;
		
		public int compareTo(Object arg0) {
			// TODO Auto-generated method stub
			if (((FieldWithDouble)arg0).v < v) return -1;
			if (((FieldWithDouble)arg0).v > v) return 1;
			return 0;
		}

	}
	

}
	