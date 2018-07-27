package Problem;

import java.util.ArrayList;

public class TwoDHWT {
	public static void ordFwdDWTForNumIters(double[][] signal, int num_iters) {
		if (signal.length != signal[0].length) {
			throw new IllegalArgumentException("Signal is not a square matrix of n x n");
		}

		if (!Utils.isPowerOf2(signal.length)) {
			throw new IllegalArgumentException("Signal's dimension is not a power of 2");
		}

		int n = signal[0].length;

		for(int i=num_iters;i>0;i--) {
			if (n > 1) {
				//Row based HWT
				applyDWTToRowsOnce(signal, n);

				//Column based HWT
				applyDWTToColsOnce(signal, n);

				//Rearranging
				fwdRearrangeAHVD(signal, n);			
			}
			n/=2;
		}
	}

	public static void applyDWTToRowsOnce(double[][] signal, int n) {
		for(int row = 0; row < n; row++) {
			for(int c = 0; c < n; c += 2) {
				double plus  =  (signal[row][c]+signal[row][c+1])/2.0;
				double minus = (signal[row][c]-signal[row][c+1])/2.0;
				signal[row][c] = plus; signal[row][c+1] = minus;
			}
		}
	}

	public static void applyDWTToColsOnce(double[][] signal, int n) {
		for(int col = 0; col < n; col++) {
			for(int r = 0; r < n; r += 2) {
				double plus  = (signal[r][col] + signal[r+1][col])/2.0;
				double minus = (signal[r][col] - signal[r+1][col])/2.0;
				signal[r][col] = plus; signal[r+1][col] = minus;
			}
		}
	}

	public static void fwdRearrangeAHVD(double[][] signal, int n) {
		double[][] M = new double[n][n];
		final int ttn_1 = n >> 1;
			int ar = 0, ac = 0;
			int hr = 0, hc = ttn_1;
			int vr = ttn_1, vc = 0;
			int dr = ttn_1, dc = ttn_1;

			for(int sig_ah_r = 0, sig_vd_r = 1; sig_ah_r < n; sig_ah_r += 2, sig_vd_r += 2) {
				ac = 0; hc = ttn_1; vc = 0; dc = ttn_1;
				for(int c = 0; c < n; c += 2) {
					M[ar][ac] = signal[sig_ah_r][c];
					M[hr][hc] = signal[sig_ah_r][c+1];
					M[vr][vc] = signal[sig_vd_r][c];
					M[dr][dc] = signal[sig_vd_r][c+1];
					ac++; hc++; vc++; dc++;
				}
				ar++; hr++; vr++; dr++;
			}
			for(int r = 0; r < n; r++) {
				System.arraycopy(M[r], 0, signal[r], 0, n);
			}
			M = null;
	}

	public static ArrayList<double[][]> orderedForwardDWTForNumIters(double[][] sample, int n, int num_iters) {
		if ( num_iters > n ) return null;
		int size = n;
		ArrayList<double[][]> iteration_rslt = null;
		ArrayList<double[][]> rslt = new ArrayList<double[][]>();
		double[][] averages = sample;
		while ( num_iters-- > 0 ) {
			iteration_rslt = orderedFastHaarWaveletTransformForNumItersAux(averages, size);
			size -= 1;
			averages = iteration_rslt.get(0);
			rslt.add(iteration_rslt.get(1));
			rslt.add(iteration_rslt.get(2));
			rslt.add(iteration_rslt.get(3));     
		}

		rslt.add(averages);
		return rslt;
	}

	public static ArrayList<double[][]> orderedFastHaarWaveletTransformForNumItersAux(double[][] sample, int n) {

		ArrayList<double[][]> rslt = new ArrayList<double[][]>();
		if ( n <= 0 ) {
			rslt.add(sample);
			rslt.add(sample);
			rslt.add(sample);
			rslt.add(sample);
			return rslt;
		}
		else {
			final int size = (int) Math.pow(2, n);
			final int half_size = size/2;
			double[][] averages = new double[half_size][half_size];
			double[][] horizontals = new double[half_size][half_size];
			double[][] verticals = new double[half_size][half_size];
			double[][] diagonals = new double[half_size][half_size];
			int i = 0, j = 0;
			double[][] tmpMatrix = new double[2][2];
			for(int tlx = 0; tlx < size; tlx += 2) {
				j = 0;
				for(int tly = 0; tly < size; tly += 2) {
					tmpMatrix[0][0] = sample[tlx][tly];
					tmpMatrix[0][1] = sample[tlx][tly+1];
					tmpMatrix[1][0] = sample[tlx+1][tly];
					tmpMatrix[1][1] = sample[tlx+1][tly+1];
					inPlaceFastHaarWaveletTransform(tmpMatrix, 2);
					averages[i][j]    = tmpMatrix[0][0];
					horizontals[i][j] = tmpMatrix[0][1];
					verticals[i][j]   = tmpMatrix[1][0];
					diagonals[i][j]   = tmpMatrix[1][1];
					j += 1;
				}
				i += 1;
			}

			rslt.add(averages);
			rslt.add(horizontals);
			rslt.add(verticals);
			rslt.add(diagonals);

			return rslt;
		}
	}

	public static void inPlaceFastHaarWaveletTransform(double[][] sample, int size) {
		inPlaceFastHaarWaveletTransformAt(sample, 0, 0, size);
	}

	public static void inPlaceFastHaarWaveletTransformAt(double[][] sample, int tlx, int tly, int quad_size) {
		if (quad_size < 2) {
			return;
		}

		if (quad_size == 2) {
			double zero_plus = sample[tlx][tly] + sample[tlx][tly + 1];
			double zero_minus = sample[tlx][tly] - sample[tlx][tly + 1];
			double one_plus = sample[tlx + 1][tly] + sample[tlx + 1][tly + 1];
			double one_minus = sample[tlx + 1][tly] - sample[tlx + 1][tly + 1];

			sample[tlx][tly] = zero_plus / 2;
			sample[tlx][tly + 1] = zero_minus / 2;
			sample[tlx + 1][tly] = one_plus / 2;
			sample[tlx + 1][tly + 1] = one_minus / 2;

			zero_plus = sample[tlx][tly] + sample[tlx + 1][tly];
			zero_minus = sample[tlx][tly] - sample[tlx + 1][tly];
			one_plus = sample[tlx][tly + 1] + sample[tlx + 1][tly + 1];
			one_minus = sample[tlx][tly + 1] - sample[tlx + 1][tly + 1];

			sample[tlx][tly] = zero_plus / 2;
			sample[tlx + 1][tly] = zero_minus / 2;
			sample[tlx][tly + 1] = one_plus / 2;
			sample[tlx + 1][tly + 1] = one_minus / 2;
		} else {
			// 1. divide the sample into four quads: Quad0, Quad1, Quad2, Quad3.
			final int half_size = quad_size / 2;
			// top left corner of Quad0
			final int quad_0_tlx = tlx;
			final int quad_0_tly = tly;
			// top left corner of Quad1
			final int quad_1_tlx = tlx;
			final int quad_1_tly = tly + half_size;
			// top left corner of Quad2
			final int quad_2_tlx = tlx + half_size;
			final int quad_2_tly = tly;
			// top left corner of Quad3
			final int quad_3_tlx = tlx + half_size;
			final int quad_3_tly = tly + half_size;

			// 2. Recursively transform each quad in place
			inPlaceFastHaarWaveletTransformAt(sample, quad_0_tlx, quad_0_tly, half_size);
			inPlaceFastHaarWaveletTransformAt(sample, quad_1_tlx, quad_1_tly, half_size);
			inPlaceFastHaarWaveletTransformAt(sample, quad_2_tlx, quad_2_tly, half_size);
			inPlaceFastHaarWaveletTransformAt(sample, quad_3_tlx, quad_3_tly, half_size);

			// 3. get the averages from each quad
			double[][] averages = new double[2][2];
			averages[0][0] = getAvrg(sample, quad_0_tlx, quad_0_tly, 0, half_size);
			averages[0][1] = getAvrg(sample, quad_0_tlx, quad_0_tly, 1, half_size);
			averages[1][0] = getAvrg(sample, quad_0_tlx, quad_0_tly, 2, half_size);
			averages[1][1] = getAvrg(sample, quad_0_tlx, quad_0_tly, 3, half_size);

			// 4. do the row transform on the averages
			sample[quad_0_tlx][quad_0_tly] = (averages[0][0] + averages[0][1]) / 2;
			sample[quad_1_tlx][quad_1_tly] = (averages[0][0] - averages[0][1]) / 2;
			sample[quad_2_tlx][quad_2_tly] = (averages[1][0] + averages[1][1]) / 2;
			sample[quad_3_tlx][quad_3_tly] = (averages[1][0] - averages[1][1]) / 2;

			// 5. do the col transform transform on the averages
			final double left_col_plus = (sample[quad_0_tlx][quad_0_tly] + sample[quad_2_tlx][quad_2_tly]) / 2;
			final double left_col_minus = (sample[quad_0_tlx][quad_0_tly] - sample[quad_2_tlx][quad_2_tly]) / 2;
			final double right_col_plus = (sample[quad_1_tlx][quad_1_tly] + sample[quad_3_tlx][quad_3_tly]) / 2;
			final double right_col_minus = (sample[quad_1_tlx][quad_1_tly] - sample[quad_3_tlx][quad_3_tly]) / 2;

			// 6. put the results back into the sample
			sample[quad_0_tlx][quad_0_tly] = left_col_plus;
			sample[quad_2_tlx][quad_2_tly] = left_col_minus;
			sample[quad_1_tlx][quad_1_tly] = right_col_plus;
			sample[quad_3_tlx][quad_3_tly] = right_col_minus;

			averages = null;
		}
	}

	public static double getAvrg(double[][] sample, int top_left_x, int top_left_y, int avrg_num, int size) {
		switch (avrg_num) {
		case 0:
			return sample[top_left_x][top_left_y];
		case 1:
			return sample[top_left_x][top_left_y + size];
		case 2:
			return sample[top_left_x + size][top_left_y];
		case 3:
			return sample[top_left_x + size][top_left_y + size];
		default:
			return -1;
		}
	}

	public static void ordInvDWTForNumIters(double[][] sig, int num_inv_iters, int num_fwd_iters) {
		if (sig.length != sig[0].length) {
			throw new IllegalArgumentException("Signal is not a square matrix of n x n");
		}

		if (!Utils.isPowerOf2(sig.length)) {
			throw new IllegalArgumentException("Signal's dimension is not a power of 2");
		}

		final int num_avail_iters = Utils.wholeLog2(sig.length);
		final int num_applied_iters = num_avail_iters - num_fwd_iters;
		int n = (int)(Math.pow(2, num_applied_iters));

		ordInvDWTForNumItersAux(sig, n, num_inv_iters, num_fwd_iters);
	}
	
	public static void ordInvDWTForNumItersAux(double[][] sig, int avrg_ttn, int num_inv_iters, int num_fwd_iters) {
        if (num_inv_iters == 0) return;
        
        invRearrangeAHVD(sig, avrg_ttn);
        for(int r = 0; r < (avrg_ttn<<1); r+=2) {
            for(int c = 0; c < (avrg_ttn<<1); c+=2) {
                invTwoByTwoDWT(sig, r, c);
            }
        }
        ordInvDWTForNumItersAuxAttic(sig, avrg_ttn<<1, num_inv_iters-1);
    }
    
    public static void ordInvDWTForNumItersAuxAttic(double[][] sig, int avrg_ttn, int num_iters) {
        if (num_iters == 0) return;
        
        invRearrangeAHVD(sig, avrg_ttn);
        for(int r = 0; r < (avrg_ttn<<1); r+=2) {
            for(int c = 0; c < (avrg_ttn<<1); c+=2) {
                invTwoByTwoDWT(sig, r, c);
            }
        }
        ordInvDWTForNumItersAuxAttic(sig, avrg_ttn<<1, num_iters-1);
    }
	
	public static void invRearrangeAHVD(double[][] sig, int avrg_ttn) {
        final int n = (avrg_ttn<<1);
        double[][] M = new double[n][n];
        for(int avrg_r=0, mr=0; avrg_r < avrg_ttn; avrg_r++, mr += 2) {
            for(int avrg_c=0, mc=0; avrg_c < avrg_ttn; avrg_c++, mc += 2) {
                M[mr][mc] = sig[avrg_r][avrg_c];
                M[mr][mc+1] = sig[avrg_r][avrg_c+avrg_ttn];
                M[mr+1][mc] = sig[avrg_r+avrg_ttn][avrg_c];
                M[mr+1][mc+1] = sig[avrg_r+avrg_ttn][avrg_c+avrg_ttn];
            }
        }
        for(int r = 0; r < n; r++) {
            System.arraycopy(M[r], 0, sig[r], 0, n);
        }
        M=null;
    }
	
	public static void invTwoByTwoDWT(double[][] sig, int top_left_r, int top_left_c) {
        final double a0 = sig[top_left_r][top_left_c]; 
        final double h0 = sig[top_left_r][top_left_c+1];
        final double v0 = sig[top_left_r+1][top_left_c]; 
        final double d0 = sig[top_left_r+1][top_left_c+1];
        final double a1_0 = a0 + h0 + v0 + d0;
        final double a1_1 = a0 - h0 + v0 - d0;
        final double a1_2 = a0 + h0 - v0 - d0;
        final double a1_3 = a0 - h0 - v0 + d0;
        sig[top_left_r][top_left_c]     = a1_0; 
        sig[top_left_r][top_left_c+1]   = a1_1;
        sig[top_left_r+1][top_left_c]   = a1_2; 
        sig[top_left_r+1][top_left_c+1] = a1_3;
    }
}