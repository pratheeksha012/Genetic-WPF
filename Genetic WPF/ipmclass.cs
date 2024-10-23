using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Policy;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
using System.Runtime.ConstrainedExecution;
using static System.Runtime.InteropServices.JavaScript.JSType;
using System.Runtime.Intrinsics;
using System.Diagnostics;
using System.Runtime.Intrinsics.Arm;
using System.Net.NetworkInformation;
using GeneticSharp;

namespace Genetic_WPF
{
    public class ipmclass
    {

        //1
        public double initial;
        public double increment;

        //4100
        public List<double> rpm;
        public List<double> we;
        public List<double> ppsi1;
        public List<double> ppsi2;
        public List<double> oLd;
        public List<double> oLq;
        public List<double> MtpvId;
        public List<double> MtpvIq;

        //1
        public double webwe;

        //4100
        public List<double> IPMstatus;
        public List<double> Id;
        public List<double> Iq;
        public List<double> Idc;
        public List<double> Irms;
        public List<double> Vd;
        public List<double> Vq;
        public List<double> Vdc;
        public List<double> Tn;
        public List<double> Pe;
        public List<double> PF;
        public List<double> P;
        public List<double> reltrq;
        public List<double> Pbtt;
        public List<double> Ploss;
        public List<double> Pcu;
        public List<double> Pfe;
        public List<double> Pstr;
        public List<double> Pfric;
        public List<double> Pwind;
        public List<double> Pinv;
        public List<double> posin;
        public List<double> negan;
        public List<double> n;
        public static List<double> Temp;
        public List<double> psipm;

        //41
        public List<double> beta;
        public List<double> Im;

        //13x4100
        public List<List<double>> tmpId;
        public List<List<double>> tmpIq;
        public List<List<double>> plaId;
        public List<List<double>> plaIq;

        //6x28  //168 flattened at some point
        public List<List<double>> Ld;
        public List<List<double>> Lq;
        public List<List<double>> psi;

        //28
        public List<double> e;
        public List<double> PsiArr;

        //6x27
        public List<List<double>> Hdmat;
        public List<List<double>> Hqmat;

        public ipmclass() {
            //initial = 0;
            //double initial;
            //double incriment;

            ////4100
            //List<double> rpm;
            //List<double> we;
            //List<double> ppsi1;
            //List<double> ppsi2;
            //List<double> oLd;
            //List<double> oLq;
            //List<double> MtpvId;
            //List<double> MtpvIq;

            //double webwe;

            ////4100
            //List<double> IPMstatus;
            //List<double> Id;
            //List<double> Iq;
            //List<double> Idc;
            //List<double> Irms;
            //List<double> Vd;
            //List<double> Vq;
            //List<double> Vdc;
            //List<double> Tn;
            //List<double> Pe;
            //List<double> PF;
            //List<double> P;
            //List<double> reltrq;
            //List<double> Pbtt;
            //List<double> Ploss;
            //List<double> Pcu;
            //List<double> Pfe;
            //List<double> Pstr;
            //List<double> Pfric;
            //List<double> Pwind;
            //List<double> Pinv;
            //List<double> posin;
            //List<double> negan;
            //List<double> n;
            //List<double> Temp;
            //List<double> psipm;

        }


        //public static Complex complex_radians_to_degrees(Complex z)
        //{

        //    double real_part = z.Real;
        //    double imaginary_part = z.Imaginary;
        //    double degrees_real = Math.Degrees(real_part);
        //    doubel degrees_imaginary = math.degrees(imaginary_part)
        //    res = complex(degrees_real, degrees_imaginary)



        //    return new Complex(real, imaginary);

        //}
        public void init(Dictionary<string, double> mot, List<List<double>> DTtab, int DTnum, Dictionary<string, List<List<List<double>>>> PLmat, double tnflag, double peflag, double _const)
        {


            List<double> temp = new List<double>();

            for (double i = initial; i < 19900.0; i += increment)
            {
                temp.Add(i);
            }

            int nrpm = temp.Count();

            List<double> prctg = new List<double>(); ;

            for (double i = 0; i <= 100.0; i += 2.5)
            {
                prctg.Add(i);

            }

            prctg[0] = 0.000001;
            int nprctg = prctg.Count();

            int num = nrpm * nprctg; //100x41 = 4100

            List<int> check = Enumerable.Repeat(1, num).ToList();

            List<int> cI = new List<int>();

            // Add the initial value (700) to the vector
            cI.Add(700);

            // Generate the range and concatenate it with the initial value
            for (int i = 650; i >= -1; i -= 25)
            {
                cI.Add(i);
            }

            List<double> wr = new List<double>(new double[num]);  // Col D
            we = new List<double>(new double[num]);  // Col E
            rpm = new List<double>(new double[num]);  // Col C
            List<double> idc = new List<double>(new double[nprctg]);
            //Idc = new List<Complex>(new Complex[num]);  // Col AU
            List<double> tXIdc = new List<double>(new double[num]);  // Col AU * Col AW

            ppsi1 = new List<double>(new double[num]);
            ppsi2 = new List<double>(new double[num]);
            oLd = new List<double>(new double[num]);
            oLq = new List<double>(new double[num]);
            List<double> dif = new List<double>(new double[num]);
            List<double> E = new List<double>(new double[num]);
            List<double> MTPAid = new List<double>(new double[num]);
            List<double> MTPAiq = new List<double>(new double[num]);
            List<double> MTPAbeta = new List<double>(new double[num]);
            List<double> MTPAIrms = new List<double>(new double[num]);
            List<double> MTPAtn = new List<double>(new double[num]);

            List<double> wea = new List<double>(new double[num]);
            List<Complex> CPSRid = new List<Complex>(new Complex[num]);
            List<Complex> CPSRiq = new List<Complex>(new Complex[num]);
            List<Complex> CPSRIrms = new List<Complex>(new Complex[num]);
            List<Complex> CPSRbeta = new List<Complex>(new Complex[num]);
            List<Complex> CPSRTn = new List<Complex>(new Complex[num]);
            List<double> XTn = new List<double>(new double[num]);
            MtpvId = new List<double>(new double[num]);
            MtpvIq = new List<double>(new double[num]);
            List<double> MtpvIm = new List<double>(new double[num]);
            List<double> MtpvIrms = new List<double>(new double[num]);
            List<double> mtpvbeta = new List<double>(new double[num]);
            List<double> MtpvTn = new List<double>(new double[num]);
            List<double> Lambda = new List<double>(new double[num]);
            List<double> Hd = new List<double>(new double[num]);
            List<double> Hq = new List<double>(new double[num]);
            List<double> H2 = new List<double>(new double[num]);
            List<double> Tlin = new List<double>(new double[num]);
            List<double> Tcal = new List<double>(new double[num]);
            List<double> kofM123 = new List<double>(new double[num]);
            List<double> kofM45 = new List<double>(new double[num]);
            List<double> kK0 = new List<double>(new double[num]);
            Tn = new List<double>(new double[num]);
            //List<double> Tn = new List<double>(new double[num]);
            List<double> Hwea = new List<double>(new double[num]);
            List<double> Hwearpm = new List<double>(new double[num]);
            //IPMstatus = new List<double>(new double[num]);
            IPMstatus = new List<double>(Enumerable.Repeat(2.0, num));
            Id = new List<double>(new double[num]);
            Iq = new List<double>(new double[num]);
            Vd = new List<double>(new double[num]);
            Vq = new List<double>(new double[num]);
            PF = new List<double>(new double[num]);
            Irms = new List<double>(new double[num]);
            P = new List<double>(new double[num]);
            reltrq = new List<double>(new double[num]);
            Pe = new List<double>(new double[num]);
            Temp = new List<double>(new double[num]);

            List<double> XVd = new List<double>(new double[num]);
            List<double> XVq = new List<double>(new double[num]);
            List<double> XVdc = new List<double>(new double[num]);
            Vdc = new List<double>(new double[num]);
            List<double> iIm = new List<double>(new double[num]);
            List<double> ibeta = new List<double>(new double[num]);
            List<double> delta = new List<double>(new double[num]);
            //List<double> Rs = new List<double>(new double[num]);
            Pcu = new List<double>(new double[num]);
            Pfe = new List<double>(new double[num]);
            Pstr = new List<double>(new double[num]);
            Pfric = new List<double>(new double[num]);
            Pwind = new List<double>(new double[num]);
            Pinv = new List<double>(new double[num]);
            Ploss = new List<double>(new double[num]);
            Pbtt = new List<double>(new double[num]);
            n = new List<double>(new double[num]);
            posin = new List<double>(new double[num]);
            negan = new List<double>(new double[num]);
            List<double> Psat = new List<double>(new double[num]);
            List<double> Pon = new List<double>(new double[num]);
            List<double> Poff = new List<double>(new double[num]);
            List<double> Pfmj = new List<double>(new double[num]);
            List<double> Prr = new List<double>(new double[num]);

            rpm = new List<double>(new double[num]);
            for (int i = 0; i < nrpm; i++)
                for (int j = 0; j < nprctg; j++)
                    rpm[i * nprctg + j] = temp[i];

            // Calculate wr (ωr)
            wr = rpm.Select(r => r / 60.0 * 2 * Math.PI).ToList(); //sus

            // Calculate we (ωe)
            we = rpm.Select(r => (r / 60.0 * 2 * Math.PI) * (mot["pole"] / 2)).ToList();

            // Calculate raw Idc array
            for (int i = 0; i < nprctg; i++)
                idc[i] = (prctg[i] * mot["cur"] / 100.0);// = PLmat["tuning"][0][0];

            // Tile Idc across RPM values
            Idc = Enumerable.Repeat(idc, nrpm).SelectMany(i => i).ToList();

            // Calculate tuned XIdc array
            List<double> tuning = new List<double>();
            for (int i = 0; i < Idc.Count(); i++)
                tuning.Add(PLmat["tuning"][i][0][0]);// = PLmat["tuning"][0][0];

            //var testString = GeneticAlgorithm.Convert1DListToString(tuning);
            //_mainWindow.myTextbox.Text = testString;

            //Console.WriteLine(Idc.Count());
            //Console.WriteLine(PLmat["tuning"][0][0].Count());

            for (int i = 0; i < Idc.Count(); i++)
                tXIdc[i] = Idc[i] * PLmat["tuning"][i][0][0];

            //for (int i = 0; i < tXIdc.Count(); i++)
            //    Console.WriteLine("{0:F10}", tXIdc[i]);

            // Initialization
            psipm = new List<double>(new double[num]);
            for (int i = 0; i < PLmat["copy"].Count(); i++)
                psipm[i] = PLmat["copy"][i][0][0];
            //psipm = new List<double>(PLmat["copy"][0][0]);
            //Console.WriteLine(PLmat["copy"].Count());

            //for (int i = 0; i < psipm.Count(); i++)
            //    Console.WriteLine(psipm[i]);

            List<double> piq1 = Enumerable.Repeat(0.0, nrpm * nprctg).ToList();//4100
            List<double> piq2 = Enumerable.Repeat(0.0, nrpm * nprctg).ToList();//4100
            List<double> psi1 = Enumerable.Repeat(0.0, nrpm * nprctg).ToList();//4100
            List<double> psi2 = Enumerable.Repeat(0.0, nrpm * nprctg).ToList();//4100

            // Initialize Ld and Lq arrays
            List<double> Ld1 = Enumerable.Repeat(mot["ld"] * mot["unit"], nrpm * nprctg).ToList();//4100
            List<double> Ld2 = Enumerable.Repeat(mot["ld"] * mot["unit"], nrpm * nprctg).ToList();//4100
            List<double> Lq1 = Enumerable.Repeat(mot["lq"] * mot["unit"], nrpm * nprctg).ToList();//4100
            List<double> Lq2 = Enumerable.Repeat(mot["lq"] * mot["unit"], nrpm * nprctg).ToList();//4100

            List<Complex> complexpsipm = new List<Complex>();
            List<Complex> complexoLd = new List<Complex>();
            List<Complex> complexE = new List<Complex>();
            List<Complex> complextXIdc = new List<Complex>();
            List<Complex> complexwe = new List<Complex>();
            List<Complex> complexdif = new List<Complex>();
            List<Complex> complexIdc = new List<Complex>();

            Complex complex1 = 1;
            Complex complex0 = 0;
            Complex complex2 = 2;
            Complex complexSQRT2 = Math.Sqrt(2);
            Complex complex180 = 180;
            Complex complexM_PI = 3.141592653589793238462643383279502884197169399375105820974944;

            int idx = 0;

            for (int iter = 0; iter < nrpm; ++iter)  // 100 iterations
            {
                //int j = 0;
                //for (int i = idx; i < idx + nprctg; ++i)
                //{
                //    ipm["rpm"][i][0][0] = temp[iter];
                //    wr[i] = (temp[iter] / 60) * 2 * M_PI;
                //    we[i] = (temp[iter] / 60 * 2 * M_PI) * (pole / 2);
                //    XIdc[i] = idc[j];
                //    tXIdc[i] = XIdc[i] * tuning[i];
                //    ++j;
                //}

                for (int ind = 0; ind < nprctg; ++ind) //41 iterations
                {

                    List<double> maxindtemp = new List<double>();
                    int maxindtempflag = 0;
                    for (int i = 0; i < cI.Count(); i++)
                    {
                        if (cI[i] <= tXIdc[ind])
                        {
                            maxindtemp.Add(i);
                            maxindtempflag = 1;
                        }
                    }
                    int? maxind = null;
                    if (maxindtempflag == 1)
                        maxind = (int)maxindtemp[0];

                    if (maxind >= 0)
                    {
                        //Console.WriteLine(idx);
                        //Console.WriteLine(maxind.Value);
                        //Console.WriteLine(" ");

                        piq1[idx] = cI[maxind.Value];            // Assign value from cI to piq1 at index idx
                        psi1[idx] = PLmat["PsiArr"][maxind.Value][0][0];// Assign value from PsiArr to ipm["ppsi1"] at index idx
                    }
                    else
                    {
                        piq1[idx] = mot["flux"];                   // Assign flux to piq1 at index idx
                        psi1[idx] = mot["flux"];          // Assign flux to ipm["ppsi1"] at index idx
                    }

                    List<double> minindtemp = new List<double>();
                    int minindtempflag = 0;
                    for (int i = 0; i < cI.Count(); i++)
                    {
                        if (cI[i] >= tXIdc[ind])
                        {
                            minindtemp.Add(i);
                            minindtempflag = 1;
                        }
                    }
                    int? minind = null;
                    if (minindtempflag == 1)
                        minind = (int)minindtemp[minindtemp.Count() - 1];

                    if (minind >= 0)
                    {
                        piq2[idx] = cI[minind.Value];            // Assign value from cI to piq2 at index idx
                        psi2[idx] = PLmat["PsiArr"][minind.Value][0][0];// Assign value from PsiArr to ipm["ppsi2"] at index idx
                    }
                    else
                    {
                        piq2[idx] = mot["flux"];                   // Assign flux to piq2 at index idx
                        psi2[idx] = mot["flux"];           // Assign flux to ipm["ppsi2"] at index idx
                    }

                    if (_const != 1)
                    {
                        if (maxind >= 0)
                        {
                            Ld1[idx] = PLmat["LdArr"][maxind.Value][0][0];
                            Lq1[idx] = PLmat["LqArr"][maxind.Value][0][0];
                        }
                        else
                        {
                            Ld1[idx] = 0;
                            Lq1[idx] = 0;
                        }

                        if (minind >= 0)
                        {
                            Ld2[idx] = PLmat["LdArr"][minind.Value][0][0];
                            Lq2[idx] = PLmat["LqArr"][minind.Value][0][0];
                        }
                        else
                        {
                            Ld2[idx] = 0;
                            Lq2[idx] = 0;
                        }
                    }

                    ppsi1[idx] = psi1[idx];
                    ppsi2[idx] = psi2[idx];
                    if (_const != 1)
                    {
                        psipm[idx] = (tXIdc[idx] - piq1[idx]) / (piq2[idx] - piq1[idx]) * (ppsi2[idx] - ppsi1[idx]) + ppsi1[idx];
                    }

                    oLd[idx] = (Idc[idx] - piq1[idx]) / (piq2[idx] - piq1[idx]) * (Ld2[idx] - Ld1[idx]) + Ld1[idx];
                    oLq[idx] = (Idc[idx] - piq1[idx]) / (piq2[idx] - piq1[idx]) * (Lq2[idx] - Lq1[idx]) + Lq1[idx];
                    dif[idx] = oLq[idx] - oLd[idx];
                    E[idx] = oLq[idx] / oLd[idx];

                    MTPAid[idx] = (psipm[idx] - Math.Sqrt(Math.Pow(psipm[idx], 2) + (8 * Math.Pow(Idc[idx], 2) * Math.Pow(dif[idx], 2)))) / (4 * dif[idx]);
                    MTPAiq[idx] = Math.Sqrt(Math.Pow(Idc[idx], 2) - Math.Pow(MTPAid[idx], 2));

                    MTPAbeta[idx] = Math.Atan2(-MTPAid[idx], MTPAiq[idx]) * 180.0 / Math.PI;                                                                          // Col Y
                    MTPAIrms[idx] = Math.Sqrt(Math.Pow(MTPAid[idx], 2) + Math.Pow(MTPAiq[idx], 2)) / Math.Sqrt(2);                                                                    // Col Z
                    MTPAtn[idx] = (3.0 / 4.0) * mot["pole"] * (psipm[idx] * MTPAIrms[idx] * Math.Cos((MTPAbeta[idx]) * Math.PI / 180.0) + dif[idx] * Math.Pow(MTPAIrms[idx], 2) * Math.Sin((2 * MTPAbeta[idx]) * Math.PI / 180.0) / 2);  // Col AA
                    wea[idx] = mot["vol"] / Math.Sqrt(Math.Pow(oLq[idx] * MTPAiq[idx], 2) + Math.Pow(oLd[idx] * MTPAid[idx] + psipm[idx], 2));  // Col AA

                    complexpsipm.Add(psipm[idx]);
                    complexoLd.Add(oLd[idx]);
                    complexE.Add(E[idx]);
                    complextXIdc.Add(tXIdc[idx]);
                    complexwe.Add(we[idx]);
                    complexdif.Add(dif[idx]);
                    complexIdc.Add(Idc[idx]);

                    CPSRid[idx] = (complexpsipm[idx] / complexoLd[idx] - Complex.Sqrt(((complexpsipm[idx] * complexpsipm[idx]) / (complexoLd[idx] * complexoLd[idx])) + ((complexE[idx] * complexE[idx]) - complex1) * ((complexpsipm[idx] * complexpsipm[idx]) / (complexoLd[idx] * complexoLd[idx]) + (complexE[idx] * complexE[idx]) * (complexIdc[idx] * complexIdc[idx]) - (mot["vol"] * mot["vol"]) / ((complexwe[idx] * complexwe[idx]) * (complexoLd[idx] * complexoLd[idx]))))) / ((complexE[idx] * complexE[idx]) - complex1);

                    //check the following
                    //complexoLd correct
                    //complexE correct
                    //complexIdc almost correct
                    //complexwe correct

                    //Console.WriteLine(CPSRid[idx]);

                    CPSRiq[idx] = Math.Sqrt(0.001);

                    if ((Complex.Pow(complextXIdc[idx], 2) - Complex.Pow(CPSRid[idx], 2)).Real > 0)
                    {

                        CPSRiq[idx] = Complex.Sqrt(Complex.Pow(complexIdc[idx], 2) - Complex.Pow(CPSRid[idx], 2));
                        //CPSRiq[idx] = sqrt((complextXIdc - CPSRid[idx])* (complextXIdc + CPSRid[idx]));
                        //CPSRiq[idx] = Complex.Sqrt((complexIdc[idx] - CPSRid[idx]) * (complexIdc[idx] + CPSRid[idx]));

                    }
                    else
                    {
                        CPSRiq[idx] = Complex.Sqrt(0.001);
                    }
                    //Console.WriteLine(CPSRiq[idx]);

                    CPSRIrms[idx] = Complex.Sqrt((CPSRid[idx] * CPSRid[idx] + CPSRiq[idx] * CPSRiq[idx]) / complex2); // Col AD

                    CPSRbeta[idx] = Complex.Atan(-CPSRid[idx] / CPSRiq[idx]) * complex180 / complexM_PI;

                    CPSRTn[idx] = (3.0 / 4.0) * mot["pole"] * (psipm[idx] * CPSRIrms[idx] * Complex.Cos((CPSRbeta[idx]) * complexM_PI / complex180) + complexdif[idx] * CPSRIrms[idx] * CPSRIrms[idx] * Complex.Sin((complex2 * CPSRbeta[idx]) * complexM_PI / complex180) / complex2);  // Col AH complex180 / complexM_PI
                    //Console.WriteLine(CPSRTn[idx]);

                    //MTPV variables
                    Lambda[idx] = ((-oLq[idx]) * psipm[idx] + Math.Sqrt((oLq[idx] * oLq[idx]) * (psipm[idx] * psipm[idx]) + (8.0 * (dif[idx] * dif[idx]) * ((mot["vol"] / we[idx]) * (mot["vol"] / we[idx]))))) / (4.0 * (-dif[idx]));  // Col AL
                    MtpvId[idx] = (Lambda[idx] - psipm[idx]) / oLd[idx]; // Col AM  #((vol / wea[idx]) * *2)))) /
                    MtpvIq[idx] = Math.Sqrt(Math.Pow((mot["vol"] / we[idx]), 2) - Math.Pow(Lambda[idx], 2)) / oLq[idx];

                    MtpvIm[idx] = Math.Sqrt(Math.Pow(MtpvId[idx], 2) + Math.Pow(MtpvIq[idx], 2));//  # Col AO
                    MtpvIrms[idx] = MtpvIm[idx] / Math.Sqrt(2);//  # Col AK
                    mtpvbeta[idx] = (Math.Atan(-MtpvId[idx] / MtpvIq[idx])) * 180.0 / Math.PI;//  # Col AJ

                    MtpvTn[idx] = (3.0 / 4.0) * mot["pole"] * (psipm[idx] * MtpvIrms[idx] * Math.Cos((mtpvbeta[idx]) * Math.PI / 180.0) + dif[idx] * MtpvIrms[idx] * MtpvIrms[idx] * Math.Sin((2 * mtpvbeta[idx]) * Math.PI / 180.0) / 2);//  # Col AP

                    check[idx] = 0;
                    if (MtpvIm[idx] < Idc[idx])
                        check[idx] = 1;

                    if (idx > 1)
                    {

                        if (check[idx] == 1)
                        {
                            if (IPMstatus[idx - 1] != 0 && IPMstatus[idx - 1] != 1)
                                IPMstatus[idx] = 1;  // MTPV
                            else
                                IPMstatus[idx] = 0;  // Skip
                        }
                        else
                            if (we[idx] < wea[idx])
                            IPMstatus[idx] = 2;  // MTPA
                        else
                            IPMstatus[idx] = 3;  // CPSR
                    }

                    Id[idx] = MtpvId[idx];
                    Iq[idx] = MtpvIq[idx];


                    if (IPMstatus[idx] == 2)
                    {
                        Id[idx] = MTPAid[idx];
                        Iq[idx] = MTPAiq[idx];
                    }
                    else if (IPMstatus[idx] == 3)
                    {
                        Id[idx] = CPSRid[idx].Real;
                        Iq[idx] = CPSRiq[idx].Real;
                    }
                    else
                    {
                        Id[idx] = MtpvId[idx];
                        Iq[idx] = MtpvIq[idx];
                    }



                    XVd[idx] = mot["res"] * Id[idx] - we[idx] * oLq[idx] * Iq[idx];  // Col AX
                    XVq[idx] = mot["res"] * Iq[idx] + we[idx] * (oLd[idx] * Id[idx] + psipm[idx]);  // Col AY
                    XVdc[idx] = Math.Sqrt(XVd[idx] * XVd[idx] + XVq[idx] * XVq[idx]);  // Col AZ

                    XTn[idx] = CPSRTn[idx].Real;  // Make a copy of CPSRTn

                    if (IPMstatus[idx] == 1)
                        XTn[idx] = MtpvTn[idx];
                    else if (we[idx] < wea[idx])
                        XTn[idx] = MTPAtn[idx];

                   // Console.WriteLine(idx + " : " + XTn[idx]);

                    //check XTn values
                    Hd[idx] = DTtab[DTnum][4] + DTtab[DTnum][2] * Id[idx] * 0.000001;  // Col BB
                    Hq[idx] = DTtab[DTnum][3] * Iq[idx] * 0.000001 + DTtab[DTnum][5];  // Col BC
                    H2[idx] = (Hd[idx] / DTtab[DTnum][0]) * (Hd[idx] / DTtab[DTnum][0]) + (Hq[idx] / DTtab[DTnum][1]) * (Hq[idx] / DTtab[DTnum][1]);  // H²: Col BD

                    Tlin[idx] = (3.0 / 2.0) * (mot["pole"] / 2.0) * (Hd[idx] * Iq[idx] - Hq[idx] * Id[idx]);  // Col BE

                    if (DTnum == 0 || DTnum == 2)
                        kofM123[idx] = DTtab[DTnum][6] * Math.Sqrt(H2[idx]) / Math.Sqrt(1 + Math.Pow(Math.Sqrt(H2[idx]) * DTtab[DTnum][6], 2));
                    else if (DTnum == 1)
                        kofM123[idx] = Math.Tanh(DTtab[DTnum][6] * Math.Sqrt(H2[idx]));
                    else
                        kofM123[idx] = 0;


                    if (DTnum == 3)
                        kofM45[idx] = Math.Atan(DTtab[DTnum][6] * Math.Sqrt(H2[idx]));
                    else if (DTnum == 4)
                        kofM45[idx] = 1 - Math.Exp(-DTtab[DTnum][6] * Math.Sqrt(H2[idx]));
                    else
                        kofM45[idx] = DTtab[DTnum][6] * Math.Sqrt(H2[idx]) / Math.Sqrt(1 + Math.Pow(Math.Sqrt(H2[idx]) * DTtab[DTnum][6], 2));


                    if (DTnum == 0 || DTnum == 1 || DTnum == 2)
                        kK0[idx] = kofM123[idx] / Math.Sqrt(H2[idx]);

                    else
                        kK0[idx] = kofM45[idx] / Math.Sqrt(H2[idx]);


                    Tcal[idx] = Tlin[idx] * kK0[idx];

                    if (IPMstatus[idx] == 2)
                        Hwea[idx] = Vdc[idx] / Math.Sqrt(Math.Pow(Hd[idx], 2) + Math.Pow(Hq[idx], 2));
                    else
                        Hwea[idx] = 0;

                    iIm[idx] = Math.Sqrt(Math.Pow(Id[idx], 2) + Math.Pow(Iq[idx], 2));
                    ibeta[idx] = Math.Atan2(-Id[idx], Iq[idx]) * 180 / Math.PI;

                    Vd[idx] = mot["res"] * Id[idx] - oLq[idx] * we[idx] * Iq[idx];
                    Vq[idx] = mot["res"] * Iq[idx] + we[idx] * (oLd[idx] * Id[idx] + psipm[idx]);
                    Vdc[idx] = Math.Sqrt(Math.Pow(Vd[idx], 2) + Math.Pow(Vq[idx], 2));

                    delta[idx] = (Math.Atan(-Vd[idx] / Vq[idx])) * 180.0 / Math.PI;

                    PF[idx] = Math.Cos((delta[idx] - ibeta[idx]) * Math.PI / 180.0);

                    Irms[idx] = iIm[idx] / Math.Sqrt(2);
                    P[idx] = (Irms[idx] * Vdc[idx]) / 1000.0;

                    if (tnflag == 1)
                        Tn[idx] = XTn[idx];

                    else
                        Tn[idx] = Tcal[idx];

                    

                    reltrq[idx] = ((3.0 / 4.0) * mot["pole"] * (dif[idx] * Math.Pow((Irms[idx] / Math.Sqrt(2)), 2) * Math.Sin(2 * ibeta[idx] * Math.PI / 180) / 2)) / Tn[idx];

                    if (peflag == 1)
                        Pe[idx] = wr[idx] * Tn[idx] / 1000.0;

                    else
                        Pe[idx] = Hwea[idx] * Tn[idx] / 1000.0;

                    //Rs[idx] = res * (1 + (ipm["Temp"][idx][0][0] - RoomT) / (234.5 + RoomT));
                    idx = idx + 1;

                }//ind
            }//iter

            //foreach (var value in Tn)
            //{
            //    Console.WriteLine(value);
            //}
            beta = ibeta.GetRange(nprctg, nprctg);
            Im = iIm.GetRange(nprctg, nprctg);

            for (int i = 0; i < nprctg; i++)
                Im[i] = Im[i] * Math.Pow(2, 0.5);


            for (int i = 0; i < check.Count; i++)
                if (check[i] == 1)
                {
                    idx = i;
                    break;
                }
            webwe = we[idx];

            Hdmat = new List<List<double>>(); //6x27
            Hqmat = new List<List<double>>(); //6x27
            Ld = new List<List<double>>(); //6x28
            Lq = new List<List<double>>(); //6x28
            psi = new List<List<double>>(); //6x28
            PsiArr = new List<double>(new double[PLmat["PsiArr"].Count]); //28
            e = new List<double>(new double[PLmat["LqArr"].Count]); //28

            // Loop to create 6 rows
            for (int i = 0; i < 6; i++)
            {
                Hdmat.Add(new List<double>(new double[27]));
                Hqmat.Add(new List<double>(new double[27]));
                Ld.Add(new List<double>(new double[28]));
                Lq.Add(new List<double>(new double[28]));
                psi.Add(new List<double>(new double[28]));
            }

            for (int i = 0; i < PLmat["Hdmat"].Count; ++i) //6
                for (int j = 0; j < PLmat["Hdmat"][0].Count; ++j) //27
                {
                    Hqmat[i][j] = PLmat["Hqmat"][i][j][0];    //6x27
                    Hdmat[i][j] = PLmat["Hdmat"][i][j][0];    //6x27
                }

            for (int i = 0; i < PLmat["Ld"].Count; ++i) //6
                for (int j = 0; j < PLmat["Ld"][0].Count; ++j) //28
                {
                    Ld[i][j] = PLmat["Ld"][i][j][0];    //6x28
                    Lq[i][j] = PLmat["Lq"][i][j][0];    //6x28
                    psi[i][j] = PLmat["Psi"][i][j][0];    //6x28
                }

            //Console.WriteLine(idx + " : " + XTn[idx]);
            //Console.WriteLine(PsiArr.Count);

            for (int i = 0; i < PLmat["PsiArr"].Count; i++) //28
            {

                PsiArr[i] = PLmat["PsiArr"][i][0][0];
                e[i] = PLmat["LqArr"][i][0][0] / PLmat["LdArr"][i][0][0];
                //Console.WriteLine(i + " : " + PsiArr[i].ToString("F8"));
                //Console.WriteLine(i + " : " + e[i].ToString("F8"));

            }

        }//init


        public void losscalc(Dictionary<string, double> mot, Dictionary<string, double> inv, Dictionary<string, double> igbt, Dictionary<string, double> Temp, Dictionary<string, double> Flag, int num)
        {

            List<double> Rs = new List<double>(new double[num]);
            List<double> iIrms = new List<double>(new double[num]);
            List<double> Mod = new List<double>(new double[num]);

            Pcu = new List<double>(new double[num]);
            Pfe = new List<double>(new double[num]);

            List<double> Psat = new List<double>(new double[num]);
            List<double> Pfmj = new List<double>(new double[num]);
            List<double> Pon = new List<double>(new double[num]);
            List<double> Poff = new List<double>(new double[num]);
            List<double> Prr = new List<double>(new double[num]);

            //Lamda expression
            //Pcu = Rs.Select((rsValue, idx) => (3.0 / 2.0 * rsValue * (Math.Pow(Id[idx], 2) + Math.Pow(Iq[idx], 2))) / 1000).ToList();

            for (int idx = 0; idx < num; idx++)
            {
                ipmclass.Temp[idx] = Temp["A"] * rpm[idx] + (Temp["B"] * Tn[idx] * Tn[idx] - Temp["C"] * Tn[idx] + Temp["iniTemp"]);
                Rs[idx] = mot["res"] * (1 + (ipmclass.Temp[idx] - Temp["RoomT"]) / (234.5 + Temp["RoomT"]));
                iIrms[idx] = Irms[idx] / Math.Sqrt(2); // Col DB
                Mod[idx] = Vdc[idx] / mot["DCV"]; // Col DC

                if (Flag["Pcu"] == 1)
                {
                    if (Flag["temp"] == 1)
                    {
                        Pcu[idx] = ((3.0 / 2.0) * Rs[idx] * (Id[idx] * Id[idx] + Iq[idx] * Iq[idx])) / 1000.0;
                    }
                    else
                    {
                        Pcu[idx] = ((3.0 / 2.0) * mot["res"] * (Id[idx] * Id[idx] + Iq[idx] * Iq[idx])) / 1000.0;
                    }
                }

                if (Flag["Pfe"] == 1)
                {
                    Pfe[idx] = (mot["cfe"] * Math.Pow(we[idx], 1.6) * ((oLd[idx] * Id[idx + 1] + psipm[idx]) * (oLd[idx] * Id[idx] + psipm[idx]) +
                        (oLq[idx] * Iq[idx]) * (oLq[idx] * Iq[idx]))) / 1000; // 1.6 from CJ28
                }
                else if (Flag["Pfe"] == 2)
                {
                    Pfe[idx] = (3.0 / (2 * 200) * (we[idx] * we[idx] * (oLd[idx] * Id[idx] + psipm[idx]) * (oLd[idx] * Id[idx] + psipm[idx]) +
                        (we[idx] * oLq[idx] * Iq[idx]) * (we[idx] * oLq[idx] * Iq[idx]))) / 1000; // 200 from WLTC E36
                }

                // Col CK
                if (Flag["Pstr"] == 1)
                {
                    Pstr[idx] = mot["cstr"] * mot["Cstrunit"] * Math.Pow(we[idx], 2) * (Math.Pow(Id[idx], 2) + Math.Pow(Iq[idx], 2)) / 1000;
                }

                // Col CL
                if (Flag["Pf"] == 1)
                {
                    Pfric[idx] = (0.000001 * mot["cPfric"] * Math.Pow(rpm[idx], 2) - 0.0003 * rpm[idx] + 1.4471) / 1000.0 * (Pe.Max() / 80.0);
                }

                // Col CM
                if (Flag["Pw"] == 1)
                {
                    Pwind[idx] = (mot["cPwind"] / 10000000) * (Math.Pow(we[idx] / (mot["pole"] / 2), 3)) / 1000;
                    //self.Pwind = (mot['cPwind'] / 10000000) * ((self.we) / (mot['pole'] / 2)) * *3 / 1000
                }

                // Col CN
                if (Flag["Pinv"] == 1)
                {
                    Pinv[idx] = 3 * ((2 * (mot["vol"]) * (Idc[idx]) * (inv["tr"] + inv["tf"]) * Math.Pow(10, -6)) +
                        (inv["von"] * (2 / Math.PI) * (Idc[idx]) * inv["ton"] * Math.Pow(10, -6)) +
                        ((Vdc[idx]) / 2 * (Idc[idx]) / 3 * inv["trr"] * Math.Pow(10, -6))) * (mot["freq"]) / 1000; // Col CS - Mode1Pinv

                }
                else if (Flag["Pinv"] == 2)
                {
                    // IGBT/FWD calculations
                    Psat[idx] = 2.0 * Math.Pow(iIrms[idx], 2) * igbt["A1"] * (1.0 / 8.0 + igbt["modulation"] * PF[idx] / (3.0 * Math.PI)) + Math.Sqrt(2) * iIrms[idx] * igbt["A0"] * (1.0 / (2.0 * Math.PI) + igbt["modulation"] * PF[idx] / 8.0); // Col CV
                    Pfmj[idx] = 2.0 * Math.Pow(iIrms[idx], 2) * igbt["B1"] * (1.0 / 8.0 - igbt["modulation"] * reltrq[idx] / (3.0 * Math.PI)) + Math.Sqrt(2) * iIrms[idx] * igbt["B0"] * (1.0 / (2.0 * Math.PI) - igbt["modulation"] * reltrq[idx] / 8.0); // Col CY
                    Pon[idx] = (igbt["C3"] * Math.Pow(Irms[idx], 3) + igbt["C2"] * Math.Pow(Irms[idx], 2) + (igbt["C1"] * Irms[idx]) + igbt["C0"]) / Math.PI * mot["freq"] * Mod[idx];// Col CW
                    Poff[idx] = (igbt["D3"] * Math.Pow(Irms[idx], 3) + igbt["D2"] * Math.Pow(Irms[idx], 2) + igbt["D1"] * Irms[idx] + igbt["D0"]) / Math.PI * mot["freq"] * Mod[idx]; // Col CX
                    Prr[idx] = (igbt["E3"] * Math.Pow(Irms[idx], 3) + igbt["E2"] * Math.Pow(Irms[idx], 2) + igbt["E1"] * Irms[idx] + igbt["E0"]) / Math.PI * mot["freq"] * Mod[idx]; // Col CZ
                    Pinv[idx] = (Psat[idx] + Pon[idx] + Poff[idx] + Pfmj[idx] + Prr[idx]) * 6 / 1000000; // Col CU - Mode2Pinv
                }

                //ipm["Ploss"][idx][0][0] = ipm["Pcu"][idx][0][0] + ipm["Pfe"][idx][0][0] + ipm["Pstr"][idx][0][0] + ipm["Pfric"][idx][0][0] + ipm["Pwind"][idx][0][0] + ipm["Pinv"][idx][0][0]; // Col CH
                //ipm["Pbtt"][idx][0][0] = ipm["Pe"][idx][0][0] + ipm["Ploss"][idx][0][0]; // Col CG
                //ipm["n"][idx][0][0] = ipm["Pe"][idx][0][0] / Pe[idx]; // Col CC
                //ipm["posin"][idx][0][0] = ipm["Pe"][idx][0][0] / ipm["Pbtt"][idx][0][0]; // Col CD
                //                                                                         // Col CE

                //ipm["negan"][idx][0][0] = (ipm["Pe"][idx][0][0] - ipm["Ploss"][idx][0][0]) / ipm["Pe"][idx][0][0];
                //if ((ipm["Pe"][idx][0][0] - ipm["Ploss"][idx][0][0]) / ipm["Pe"][idx][0][0] <= 0)
                //    ipm["negan"][idx][0][0] = 0;


            }
        }

    }//class
}//namespace
