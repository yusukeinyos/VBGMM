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
            double[] alpha, beta, gamma = new double[M];
            double[][] u;
            double[][,] V;
            init();
            //CsvFileIO.CsvFileIO.WriteData("motoData.csv", X);

            Update(out alpha, out beta, gamma, out u, out V);
            predicted_Distribution_output(u, V);
        }
        //----------------------------------------------------------------------------------
        static void init()
        {
            X = CsvFileIO.CsvFileIO.ReadData("GaussianMixtureData20130501.txt");
            S = X.GetLength(0);
            D = X.GetLength(1);
            M = 4;
            alpha0 = 3.0;
            beta0 = 3.0;
            gamma0 = 3.0;
            u0 = new double[] { 1.0, 1.0 };
            V0 = Mt.Mul(0.001, Mt.UnitMatrix(D));
            V0_inv = V0.Inverse();
        }
        //----------------------------------------------------------------------------------
        static void predicted_Distribution_output(double[][] u, double[][,] V)
        {
            int Num = 300; //等高線描画用プロット数
            int Ratio = 10000000;
            RandomMT rand = new RandomMT();

            double[][] x_toukou = new double[Num][]; //等高線描画用プロット点

            for (int m = 0; m < M; m++)
            {
                double[,] U; //直交固有ベクトル行列
                double[] lambda; //固有値

                Tuple<double[,], double[]> tu;
                tu = Mt.EigenValueDecompositionM(V[m]);
                lambda = tu.Item2;
                U = tu.Item1;

                bool sign_flag = true;
                for (int n = 0; n < Num; n++)
                {
                    sign_flag = !sign_flag;
                    double[] y = new double[2];
                    y[1] = 2.0 * Math.Sqrt(lambda[0] * Ratio) * rand.Double() - Math.Sqrt(lambda[0] * Ratio); //楕円の定義域内に
                    if (sign_flag)
                        y[0] = Math.Sqrt(lambda[1] * Ratio - y[1] * y[1] * lambda[1] / lambda[0]); //次元を逆にしたら上手く等高線が書けたが．．．たぶん直交行列Uの構成の仕方が原因
                    else
                        y[0] = -Math.Sqrt(lambda[1] * Ratio - y[1] * y[1] * lambda[1] / lambda[0]);


                    x_toukou[n] = Mt.Add(Mt.Mul(U.Inverse(), y), u[m]);
                }

                CsvFileIO.CsvFileIO.WriteData("x_toukou(M=" + m + ").csv", jagTomatrix(x_toukou));
            }

            CsvFileIO.CsvFileIO.WriteData("mu_matrix.csv", jagTomatrix(u));

        }
        //----------------------------------------------------------------------------------
        static void Update(out double[] alpha, out double[] beta, double[] gamma, out double[][] u, out double[][,] V)
        {
            RandomMT rnd = new RandomMT();
            double[,] lo = new double[S, M];
            double[,] r = new double[S, M];
            double[] ita = new double[M];
            alpha = new double[M]; beta = new double[M]; gamma = new double[M];
            u = new double[M][];
            V = new double[M][,];
            for (int m = 0; m < M; m++)
            {
                alpha[m] = alpha0;
                beta[m] = beta0;
                gamma[m] = gamma0;
                //u[m] = (double[])u0.Clone();
                V[m] = (double[,])V0.Clone();
            }
            u = New.Array(M, m => matrixTojag(X)[rnd.Int(S)].ToArray());

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
                    {
                        if (lo_sum != 0)
                            r[i, j] = lo[i, j] / lo_sum;
                    }
                }
                //Update ita
                for (int j = 0; j < M; j++)
                    ita[j] = sumOfmatrix_column(r, j);
                //Update alpha,beta,gamma
                old_alpha = (double[])alpha.Clone();
                old_beta = (double[])beta.Clone();
                old_gamma = (double[])gamma.Clone();


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
        static double[,] jagTomatrix(double[][] input)
        {
            double[,] output = new double[input.GetLength(0), input[0].Length];
            for (int i = 0; i < input.GetLength(0); i++)
            {
                for (int j = 0; j < input[0].Length; j++)
                {
                    output[i, j] = input[i][j];
                }
            }
            return output;
        }
        static double[][] matrixTojag(double[,] input)
        {
            return New.Array(input.GetLength(0), m => New.Array(input.GetLength(1), s => input[m, s]));
        }
    }
}
