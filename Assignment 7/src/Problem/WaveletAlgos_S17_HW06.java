package Problem;

import java.util.ArrayList;

import org.opencv.core.CvType;
import org.opencv.core.Mat;
import org.opencv.highgui.Highgui;
import org.opencv.imgproc.Imgproc;

public class WaveletAlgos_S17_HW06 {

	// Change this variable accordingly
	static final String OPENCV_DLL_PATH= "D:\\workspace\\USU-Assignments\\CS 6810 - Wavelets and Wavelet Algorithms\\Assignment 7\\Assignment 7\\src\\Problem\\opencv_java2413.dll";

	// Use a static code block to load the dll/so;
	static {
		System.load(OPENCV_DLL_PATH);
	}

	static double[][] signal_2x2_1 = {
			{9, 7},
			{5, 3}
	};

	static double[][] signal_2x2_2 = {
			{1, 0},
			{0, 1}
	};


	static double[][] signal_2x2_3 = {
			{255, 100},
			{50, 250}
	};

	static double[][] signal_2x2_4 = {
			{7, 5},
			{3, 1}
	};

	static double[][] signal_4x4_1 = {
			{9,	7,	6,	2},
			{5,	3,	4,	4},
			{8,	2,	4,	0},
			{6,	0,	2,	2}
	};

	static double[][] signal_4x4_2 = {
			{6.0,	8.0,	10.0,	8.0},
			{10.0,	4.0,	8.0,	2.0},
			{2.0,	6.0,	0.0,	6.0},
			{10.0,	6.0,	6.0,	10.0}
	};

	static double[][] signal_8x8_1 = {
			{255, 0, 0, 0, 0, 0, 0, 100},
			{0, 255, 0, 0, 0, 0, 100, 0},
			{0, 0, 255, 0, 0, 100, 0, 0},
			{0, 0, 0, 250, 100, 0, 0, 0},
			{0, 0, 0, 120, 150, 0, 0, 0},
			{0, 0, 120, 0, 0, 150, 0, 0},
			{0, 120, 0, 0, 0, 0, 150, 0},
			{120, 0, 0, 0, 0, 0, 0, 150}
	};

	static double[][] signal_8x8_2 = {
			{8,  5,  4,  8,  6, 8, 10,  8},
			{8, 10, 10,  4, 10, 4,  8,  2},
			{6, 10,  2,  4,  2, 6,  1,  6},
			{0,  8,  0, 10, 10, 6,  6, 10},
			{8,  8,  8,  0,  4, 8,  6,  2},
			{10, 6,  2,  2,  6, 6,  6,  8},
			{2,  4, 10, 10, 10, 4,  6, 10},
			{10, 6, 10,  6,  6, 4,  4,  4}
	};

	static double[][] signal_8x8_3 = {
			{255, 0, 0, 0, 0, 0, 0, 0},
			{0, 255, 0, 0, 0, 0, 0, 0},
			{0, 0, 255, 0, 0, 0, 0, 0},
			{0, 0, 0, 255, 0, 0, 0, 0},
			{0, 0, 0, 0, 255, 0, 0, 0},
			{0, 0, 0, 0, 0, 255, 0, 0},
			{0, 0, 0, 0, 0, 0, 255, 0},
			{0, 0, 0, 0, 0, 0, 0, 255}
	};

	public static double[][] img_sample_32x32_01 = {
			{255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255},
			{0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0},
			{0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0},
			{0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0},
			{0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0},
			{0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0},
			{0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0},
			{0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0},
			{255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255},
	};


	public static void createScaledImageOf2DHWT(String infile, String outfile, int num_iters) {
		double[][] mat = getGrayscalePixMat(infile);
		int n = (int)(Math.log(mat.length)/Math.log(2));

		ArrayList<double[][]> transform = TwoDHWT.orderedForwardDWTForNumIters(mat, n, num_iters);

		double[][] lastAvrgMat = transform.remove(transform.size()-1);
		Mat grayscale = new Mat(mat.length, mat[0].length, CvType.CV_8UC1);

		int avrg_size = mat.length;
		for(int i = 0, mat_size = mat.length; i < transform.size(); i += 3, mat_size /= 2) {
			scaleSetHorMatValsInImage(grayscale, mat_size, transform.get(i));
			scaleVerMatValsInImage(grayscale, mat_size, transform.get(i+1));
			scaleDigMatValsInImage(grayscale, mat_size, transform.get(i+2));
			avrg_size /= 2;
		}

		setAvrgMatValsInImage(grayscale, 2*avrg_size, lastAvrgMat);

		Highgui.imwrite(outfile, grayscale);
		grayscale.release();
	}    

	public static void testOrd2DHWT(double[][] data, int num_iters) {
		final int dim = data.length;
		TwoDHWT.ordFwdDWTForNumIters(data, num_iters);
		System.out.println("Result Matrix");
		Utils.display2DArray(data, dim, dim);
		System.out.println();
	}

	public static void main(String[] args) {
		//For 2D HWT
		test_ordFwdInvHWT2(signal_4x4_1, 2, 2);
	}

	public static double[][] getGrayscalePixMat(String infile) {
		Mat orig = Highgui.imread(infile);
		if (orig.rows() == 0 || orig.cols() == 0) {
			throw new IllegalArgumentException("Failed to read " + infile);
		}

		Mat grayscale = new Mat(orig.rows(), orig.cols(), CvType.CV_8UC1);
		Imgproc.cvtColor(orig, grayscale, Imgproc.COLOR_RGB2GRAY);

		double[][] pix_mat = get1CPixMat(grayscale);
		orig.release();
		grayscale.release();
		return pix_mat;
	}

	public static double[][] get1CPixMat(Mat orig) {
		if (orig.rows() == 0 || orig.cols() == 0) {
			throw new IllegalArgumentException("empty image");
		}

		double[][] pix_mat = new double[orig.rows()][orig.cols()];
		for (int row = 0; row < orig.rows(); row++) {
			for (int col = 0; col < orig.cols(); col++) {
				pix_mat[row][col] = orig.get(row, col)[0];
			}
		}
		return pix_mat;
	}

	public static void scaleSetHorMatValsInImage(Mat grayscale, int mat_size, double[][] cmat) {
		double max    = findMaxVal(cmat);
		double min    = findMinVal(cmat);
		double crange = max - min;
		double scaler = 255/crange;
		int half = mat_size/2;
		Mat temp = new Mat(cmat.length, cmat[0].length, CvType.CV_8UC1);
		for(int img_row = 0, cmat_row = 0; img_row < half; img_row++, cmat_row++) {
			for(int img_col = half, cmat_col = 0; img_col < mat_size; img_col++, cmat_col++) {
				double cmat_val = (cmat[cmat_row][cmat_col] - min)*scaler;
				double[] data = { cmat_val, cmat_val, cmat_val };
				grayscale.put(img_row, img_col, data);
				temp.put(cmat_row, cmat_col, data);
			}
		}
		temp.release();

	}

	public static double findMaxVal(double[][] cmat) {
		double max = Double.MIN_VALUE;
		for(int row = 0; row < cmat.length; row++) {
			for(int col = 0; col < cmat[row].length; col++) {
				if ( cmat[row][col] > max ) {
					max = cmat[row][col];
				}
			}
		}
		return max;
	}

	public static double findMinVal(double[][] cmat) {
		double min = Double.MAX_VALUE;
		for(int row = 0; row < cmat.length; row++) {
			for(int col = 0; col < cmat[row].length; col++) {
				if ( cmat[row][col] < min ) {
					min = cmat[row][col];
				}
			}
		}
		return min;
	}

	public static void scaleVerMatValsInImage(Mat grayscale, int mat_size, double[][] cmat) {
		double max    = findMaxVal(cmat);
		double min    = findMinVal(cmat);
		double crange = max - min;
		double scaler = 255/crange;
		int half = mat_size/2;
		Mat temp = new Mat(cmat.length, cmat[0].length, CvType.CV_8UC1);
		for(int img_row = half, cmat_row = 0; img_row < mat_size; img_row++, cmat_row++) {
			for(int img_col = 0, cmat_col = 0; img_col < half; img_col++, cmat_col++) {
				double cmat_val = (cmat[cmat_row][cmat_col] - min)*scaler;
				double[] data = { cmat_val, cmat_val, cmat_val };
				grayscale.put(img_row, img_col, data);
				temp.put(cmat_row, cmat_col, data);
			}
		}
		temp.release();
	}

	public static void scaleDigMatValsInImage(Mat grayscale, int mat_size, double[][] cmat) {
		double max    = findMaxVal(cmat);
		double min    = findMinVal(cmat);
		double crange = max - min;
		double scaler = 255/crange;
		int half = mat_size/2;
		Mat temp = new Mat(cmat.length, cmat[0].length, CvType.CV_8UC1);
		for(int img_row = half, cmat_row = 0; img_row < mat_size; img_row++, cmat_row++) {
			for(int img_col = half, cmat_col = 0; img_col < mat_size; img_col++, cmat_col++) {
				double cmat_val = (cmat[cmat_row][cmat_col] - min)*scaler;
				double[] data = { cmat_val, cmat_val, cmat_val };
				grayscale.put(img_row, img_col, data);
				temp.put(cmat_row, cmat_col, data);
			}
		}
		temp.release();
	}

	public static void setAvrgMatValsInImage(Mat grayscale, int mat_size, double[][] cmat) {
		double max    = findMaxVal(cmat);
		double min    = findMinVal(cmat);
		double crange = max - min;
		double scaler = 255/crange;
		int half = mat_size/2;
		Mat temp = new Mat(cmat.length, cmat[0].length, CvType.CV_8UC1);
		for(int img_row = 0, cmat_row = 0; img_row < half; img_row++, cmat_row++) {
			for(int img_col = 0, cmat_col = 0; img_col < half; img_col++, cmat_col++) {
				//double cmat_val = (cmat[cmat_row][cmat_col] - min)*scaler;
				double cmat_val = cmat[cmat_row][cmat_col];
				double[] data = { cmat_val, cmat_val, cmat_val };
				grayscale.put(img_row, img_col, data);
				temp.put(cmat_row, cmat_col, data);
			}
		}
		temp.release();
	}

	static void test_ordFwdInvHWT2(double[][] sig, int num_fwd_iters, int num_inv_iters) {
		TwoDHWT.ordFwdDWTForNumIters(sig, num_fwd_iters);
		System.out.println("================");
		System.out.println("Forward HWT Transform");
		Utils.display2DArray(sig, sig.length, sig[0].length);
		
		TwoDHWT.ordInvDWTForNumIters(sig, num_inv_iters, num_fwd_iters);
		System.out.println("================");
		System.out.println("Inverse HWT Transform");
		Utils.display2DArray(sig, sig.length, sig[0].length);
	}
}