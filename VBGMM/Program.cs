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
        static int M; //混合数
        static double alpha0;
        static double beta0;
        static double gamma0;
        static double[] u0;
        static double[,] V0;
        static double[,] V0_inv;

        static double[,] X; //サンプルデータ　(S×D)
        static void Main(string[] args)
        {
        }
        //----------------------------------------------------------------------------------
        static void init()
        {
            S = X.GetLength(0);
            D = X.GetLength(1);
        }
        //----------------------------------------------------------------------------------
        static void Update(double[] alpha, double[] beta, double[] gamma, double[][] u, double[][,] V)
        {
            double[,] lo = new double[S, M];
            double[,] r = new double[S, M];
            double[] ita = new double[M];
            //Update lo
            for (int s = 0; s < S; s++)
                for (int m = 0; m < M; m++)
                {

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
            double[,] A = new double[D, D];
            double[,] B = new double[D, D];
            double[,] rX_adamar = new double[S, D];
            for (int j = 0; j < M; j++)
            {
                for (int i = 0; i < D; i++)
                    for (int ii = 0; ii < D; ii++)
                        A[i, ii] = beta0 * (u[j][i] - u0[i]) * (u[j][ii] - u0[ii]);
                for (int s = 0; s < S; s++)
                    for (int d = 0; d < D; d++)
                        rX_adamar[s, d] = X[s, d] * r[s, j];
                B = Mt.Mul(X.T(), rX_adamar);
                V[j] = Mt.Add(Mt.Add(V0_inv, A), B);
                V[j] = V[j].Inverse();
            }

        }
        //----------------------------------------------------------------------------------
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
