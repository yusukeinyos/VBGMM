using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using InoueLab;
using CsvFileIO;

namespace VBGMM
{
    class Program
    {
        static int S; //サンプル数
        static int D; //次元数
        static int M; //最大混合数
        static double alpha0;
        static double beta0;
        static double gamma0;
        static double[] u0;
        static double[,] V0;
        static double[,] V0_inv;

        static double[,] X; //サンプルデータ　(S×D)
        static void Main(string[] args)
        {
            double[] alpha, beta, gamma=new double[M];
            double[][] u;
            double[][,] V;
            Update(out alpha, out beta, gamma, out u, out V);
        }
        //----------------------------------------------------------------------------------
        static void init()
        {
            S = X.GetLength(0);
            D = X.GetLength(1);
            M = 5;
            alpha0 = 3.0;
            beta0 = 3.0;
            gamma0 = 3.0;

        }
        //----------------------------------------------------------------------------------
        static void Update(out double[] alpha, out double[] beta, double[] gamma, out double[][] u, out double[][,] V)
        {
            double[,] lo = new double[S, M];
            double[,] r = new double[S, M];
            double[] ita = new double[M];
            alpha = new double[M]; beta = new double[M]; gamma = new double[M];
            u = new double[M][];
            V = new double[M][,];
            double[] old_alpha = new double[M], old_beta = new double[M], old_gamma = new double[M];
            var old_V = (double[][,])V.Clone();

            double error;
            do
            {
                //Update lo
                for (int s = 0; s < S; s++)
                    for (int m = 0; m < M; m++)
                    {
                        double a = diGamma(alpha[m]);
                        double b = Enumerable.Range(1, D).Sum(dd => diGamma((gamma[m] + 1 - dd) / 2.0));
                        double c = Math.Log(Mt.Det(V[m]));
                        double d = -D / beta[m];
                        double e = -gamma[m] * Mt.QuadraticForm(V[m], Mt.Sub(getRow(X, s), u[m]));
                        lo[s, m] = Math.Exp(a + 0.5 * (b + c + d + e));
                    }
                //Update r
                for (int i = 0; i < S; i++)
                {
                    double lo_sum = sumOfmatrix_row(lo, i);
                    for (int j = 0; j < M; j++)
                        r[i, j] = lo[i, j] / lo_sum;
                }
                //Update ita
                for (int j = 0; j < M; j++)
                    ita[j] = sumOfmatrix_column(r, j);
                //Update alpha,beta,gamma
                alpha.CopyTo(old_alpha);
                beta.CopyTo(old_beta);
                gamma.CopyTo(old_gamma);
                for (int j = 0; j < M; j++)
                {
                    alpha[j] = alpha0 + ita[j];
                    beta[j] = beta0 + ita[j];
                    gamma[j] = gamma0 + ita[j];
                }
                //Update u
                for (int j = 0; j < M; j++)
                {
                    double[,] rX = Mt.Mul(X.T(), r);
                    for (int i = 0; i < D; i++)
                        u[j][i] = (beta0 * u0[i] + rX[i, j]) / beta[j];
                }
                //Update V
                V.CopyTo(old_V);
                double[,] A = new double[D, D];
                double[,] B = new double[D, D];
                double[,] rXsubu_adamar = new double[S, D];
                double[,] Xsubu = new double[S, D];
                for (int j = 0; j < M; j++)
                {
                    for (int i = 0; i < D; i++)
                        for (int ii = 0; ii < D; ii++)
                            A[i, ii] = beta0 * (u[j][i] - u0[i]) * (u[j][ii] - u0[ii]);
                    for (int s = 0; s < S; s++)
                        for (int d = 0; d < D; d++)
                        {
                            rXsubu_adamar[s, d] = (X[s, d] - u[j][d]) * r[s, j];
                            Xsubu[s, d] = X[s, d] - u[j][d];
                        }
                    B = Mt.Mul(Xsubu.T(), rXsubu_adamar);
                    V[j] = Mt.Add(Mt.Add(V0_inv, A), B);
                    V[j] = V[j].Inverse();
                }
                error = errorCalc(alpha, old_alpha) + errorCalc(beta, old_beta) + errorCalc(gamma, old_gamma) + errorCalc(V, old_V);
            } while (error > 1e-10);

        }
        //----------------------------------------------------------------------------------
        static int select_mixnum()
        {
            double gamma_alpha0 = Mt.Gamma(alpha0);
            double gamma_halfgamma0 = Mt.Gamma(beta0 / 2.0);
            double det_V0 = V0.Det();

            int M_estimated = M;
            for (int m = M; m > 0; m--)
            {

            }
            return M_estimated;
        }
        //----------------------------------------------------------------------------------
        //ディガンマ関数
        static double diGamma(double x)
        {
            double delta = 1e-6;
            return (Mt.LogGamma(x + delta) - Mt.LogGamma(x)) / delta;
        }
        //----------------------------------------------------------------------------------
        static double errorCalc(double[] a, double[] b)
        {
            return Mt.InnerProduct(Mt.Sub(a, b), Mt.Sub(a, b)) / a.Length;
        }
        static double errorCalc(double[][] a, double[][] b)
        {
            return New.Array(a.Length, i => errorCalc(a[i], b[i])).Sum();
        }
        static double errorCalc(double[,] a, double[,] b)
        {
            return Mt.Det(Mt.Sub(a, b)) * Mt.Det(Mt.Sub(a, b)) / (a.GetLength(0) * a.GetLength(1));
        }
        static double errorCalc(double[][,] a, double[][,] b)
        {
            return New.Array(a.Length, i => errorCalc(a[i], b[i])).Sum();
        }
        static double[] getRow(double[,] matrix, int row_num)
        {
            double[] output = new double[matrix.GetLength(1)];
            for (int i = 0; i < matrix.GetLength(1); i++)
                output[i] = matrix[row_num, i];
            return output;
        }
        static double sumOfmatrix_row(double[,] mat, int row_num)
        {
            double sum = 0;
            for (int i = 0; i < mat.GetLength(1); i++)
                sum += mat[row_num, i];
            return sum;
        }
        static double sumOfmatrix_column(double[,] mat, int column_num)
        {
            double sum = 0;
            for (int i = 0; i < mat.GetLength(0); i++)
                sum += mat[i, column_num];
            return sum;
        }
    }
}
