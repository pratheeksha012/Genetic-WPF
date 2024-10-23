using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Intrinsics.Arm;
using System.Text;
using System.Threading.Tasks;

namespace Genetic_WPF
{



    internal class evsimclass
    {
        public double[] rpm { get; set; }   //1x650
        public double[,] T { get; set; }  //27x91
        public double[] maxT { get; set; }   //1x27
        public double[] b { get; set; }    //1x27    
        public double[,] ptgId { get; set; }  //4x91
        public double[,] ptgIq { get; set; }  //4x91
        public double? MTPATe { get; set; }  //double
        public double? MTPAwe { get; set; } //double // Variable to store MTPA ωe value
        public double? MTPAImax { get; set; }   //double
        public double? MTPArpm { get; set; }   //double
        public double? MaxSpdRpm { get; set; } //double   // Max speed rpm value obtained from EVtab
        public double[] status { get; set; } //1x650 // EVtab status variable
        public double[] absearch { get; set; } //1x650 // EVtab absearch variable
        public double[] trqTem { get; set; } // 1x13 // Torque array that will store all the 13 variations of Tem
        public double[] mrpm { get; set; } //1x3 // Combined array of maxspeed, const torque and mtpv rpm
        public double[] vLimId { get; set; }  //1x61
        public double[] vLimIq { get; set; }  //1x61
        public double[] vLimTn { get; set; }  //1x61  
        public double[] vLim_iq { get; set; }  //1x67
        public double[] vLim_I { get; set; }  //1x67
        public double[] MaxPwrId { get; set; }  //1x61
        public double[] MaxPwrIq { get; set; }   //1x61
        public double[] MaxPwrTn { get; set; }   //1x61
        public double[] Maxpwr_iq { get; set; }  //1x67
        public double[] MaxPwr_I { get; set; }   //1x67
        public double[] MTPVId { get; set; } //1x61
        public double[] MTPVIq { get; set; } //1x61
        public double[] MTPVTn { get; set; } //1x61
        public double[] MTPV_iq { get; set; }    //1x67
        public double[] MTPV_I { get; set; } //1x67
        public double[] PfId { get; set; }   //1x67
        public double[] PfIq { get; set; }   //1x67

        public evsimclass()
        {
            rpm = new double[650];
            T = new double[27, 91];
            maxT = new double[27];
            b = new double[27];
            ptgId = new double[4, 91];
            ptgIq = new double[4, 91];
            MTPATe = null;
            MTPAwe = null; // Variable to store MTPA ωe value
            MTPAImax = null;
            MTPArpm = null;
            MaxSpdRpm = null; // Max speed rpm value obtained from EVtab
            status = new double[650]; // EVtab status variable
            absearch = new double[650]; // EVtab absearch variable
            trqTem = new double[13]; // Torque array that will store all the 13 variations of Tem
            mrpm = new double[3]; // Combined array of maxspeed, const torque, and mtpv rpm
            vLimId = new double[61];
            vLimIq = new double[61];
            vLimTn = new double[61];
            vLim_iq = new double[67];
            vLim_I = new double[67];
            MaxPwrId = new double[61];
            MaxPwrIq = new double[61];
            MaxPwrTn = new double[61];
            Maxpwr_iq = new double[67];
            MaxPwr_I = new double[67];
            MTPVId = new double[61];
            MTPVIq = new double[61];
            MTPVTn = new double[61];
            MTPV_iq = new double[67];
            MTPV_I = new double[67];
            PfId = new double[67];
            PfIq = new double[67];
        }


        public void init()
        {

            // Equivalent of np.arange(0, 91)
            int[] deg = new int[91];
            for (int i = 0; i < 91; i++)
                deg[i] = i;

            // Equivalent of np.arange(0, 651, 25)
            int[] ind = new int[27];  // 0, 25, ..., 650
            for (int i = 0; i < 27; i++)
                ind[i] = i * 25;

            ind[0] = 5;
            ind[ind.Length - 1] = 700;

            // Calculate percentage (ptg = ind / ind[-1])
            double[] ptg = new double[ind.Length];
            for (int i = 0; i < ind.Length; i++)
                ptg[i] = ind[i] / (double)ind[ind.Length - 1];

            int numi = ind.Length;
            int numd = deg.Length;

            // Initialize arrays T, maxT, and b
            T = new double[numi, numd];
            maxT = new double[numi];
            b = new double[numi];

            // Calculate torque and find max torque and corresponding angle
        //    for (int idx = 0; idx < numi; idx++)
        //    {
        //        for (int dind = 0; dind < numd; dind++)
        //        {
        //            double pole = mot["pole"];
        //            double flux = mot["flux"];
        //            double DCC = mot["DCC"];
        //            double d = mot["d"];
        //            double ptgVal = ptg[idx];
        //            double degRad = Math.PI * deg[dind] / 180.0; // Convert degrees to radians

        //            T[idx, dind] = 3.0 / 4.0 * pole * (
        //                (flux * (ptgVal * DCC) * Math.Cos(degRad)) +
        //                (d * Math.Pow((ptgVal * DCC), 2) * Math.Sin(2 * degRad)) / 2.0
        //            );
        //        }
        //        MaxT[idx] = MaxArray(T, idx); // Get maximum value from T[idx, :]
        //        B[idx] = deg[ArgMaxArray(T, idx)]; // Get index of maximum value and use it to get degree
        //    }

        //    // Calculate Id, Iq at 25%, 50%, 75%, 100%
        //    double[] arr = new double[4];
        //    for (int i = 0; i < 4; i++)
        //        arr[i] = ((25 * (i + 1)) / 100.0) * 480 * Math.Sqrt(2);

        //    PtgId = new double[4, numd];
        //    PtgIq = new double[4, numd];

        //    for (int i = 0; i < 4; i++)
        //    {
        //        for (int j = 0; j < numd; j++)
        //        {
        //            double degRad = Math.PI * deg[j] / 180.0;
        //            PtgId[i, j] = -arr[i] * Math.Sin(degRad);
        //            PtgIq[i, j] = arr[i] * Math.Cos(degRad);
        //        }
        //    }
        //}

        //// Function to get the max value of the 1D array T[idx,:]
        //private double MaxArray(double[,] array, int row)
        //{
        //    double maxVal = array[row, 0];
        //    for (int i = 1; i < array.GetLength(1); i++)
        //    {
        //        if (array[row, i] > maxVal)
        //            maxVal = array[row, i];
        //    }
        //    return maxVal;
        //}

        //// Function to get the index of the max value in T[idx,:]
        //private int ArgMaxArray(double[,] array, int row)
        //{
        //    int maxIdx = 0;
        //    double maxVal = array[row, 0];
        //    for (int i = 1; i < array.GetLength(1); i++)
        //    {
        //        if (array[row, i] > maxVal)
        //        {
        //            maxVal = array[row, i];
        //            maxIdx = i;
        //        }
        //    }
        //    return maxIdx;

        //    //    }

        }
    }
    }
