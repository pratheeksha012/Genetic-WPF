using System;
using System.Collections.Generic;
using System.Drawing.Printing;
using System.Linq;
using System.Reflection.Metadata;
using System.Runtime.Intrinsics.Arm;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Threading;

namespace Genetic_WPF
{
    public class wltcclass
    {
        //1
        public double WLTCconpow;
        public double WLTCdist;
        public double bttcon;
        public double crusdist;
        public double Fgravity;
        public double num;
        public double count;
        //1801
        public List<double> sec;
        public List<double> acckm;
        public List<double> wwrpm;
        public List<double> Frolling;
        public List<double> FAero;
        public List<double> Fx;
        public List<double> werpm;
        public List<double> Tn;
        public List<double> condPe;
        public List<double> AccPe;
        public List<double> AccPbatt;
        public List<double> Temp;

        public List<double> Id;
        public List<double> Iq;
        public List<double> Pekwh;
        public List<double> Pekws;
        public List<double> Ploss;
        public List<double> Pcu;
        public List<double> Pfe;
        public List<double> Pstr;
        public List<double> Pf;
        public List<double> Pw;

        public List<double> Pinv;
        public List<double> n;
        public List<double> cdcI;
        // public List<double> Fgravity;
        //self.Temp[0] = 57.92 

        //4100
        public List<double> rpmflag;


        //ctId,ctIq

        public void init(double pole, Dictionary<string, double> btt, Dictionary<string, double> gear, Dictionary<string, double> tire, Dictionary<string, dynamic>ip, ipmclass ipm)
        {
            double angle = (Math.Atan(ip["rdgrd"]) * 180 / Math.PI);
            double slopesin = Math.Sin(angle * (Math.PI / 180.0));   
            double slopecos = Math.Cos(angle * (Math.PI / 180.0));

            double Fgravity = (ip["GrossWt"]) * (ip["gravity"]) * slopesin;

            //4100
            List<double> rpmflag= new List<double>();
            List<double> tmprpm = new List<double>();
            List<bool> idx1 = new List<bool>();

            ipm.rpm.Add(0);
            ipm.Idc.Add(0);
            for (int i = 0; i < 4100; i++)
            {
                //tmprpm.Add(ipm.rpm[i] - ipm.increment);
                if ((ipm.rpm[i] < ipm.rpm[i + 1]) && (ipm.Idc[i + 1] < ipm.Idc[i]))
                {
                    rpmflag.Add(0.0);
                }
                else
                {
                    rpmflag.Add(1.0);
                }

                // Initialize with default value (0.0)
                tmprpm.Add(ipm.rpm[i] - ipm.increment);
                if (i == 0)
                {
                    tmprpm.Add(0);
                }
                //
                //  if (i < 4098)
                //  {
                //      idx1.Add((ipm.rpm[i] < ipm.rpm[i + 1]) && (ipm.Idc[i + 1] < ipm.Idc[i]));
                //  }
                //  else
                //   {
                //      idx1.Add(false);
                //   }

                //   if (idx1[i] == true)
                //   {
                //       rpmflag.Add(0);
                //   }
                // Console.WriteLine(i);
                //Console.WriteLine(rpmflag[i]);
            }
            ipm.rpm.RemoveAt(ipm.rpm.Count - 1);
            ipm.Idc.RemoveAt(ipm.Idc.Count - 1);
            //1801
            List<double> acckm = new List<double>();
            List<double> InertialL = new List<double>();
            List<double> acc = new List<double>();
            List<double> dist = new List<double>();
            List<double> wwrpm = new List<double>();
            List<double> ww = new List<double>();
            List<double> alpha = new List<double>();
            List<double> Frolling = new List<double>();
            List<double> FAero = new List<double>();
            List<double> Fx = new List<double>();
            List<double> Tw = new List<double>();
            List<double> Pe = new List<double>();
            List<double> we = new List<double>();
            List<double> wr = new List<double>();
            List<double> werpm = new List<double>();
            List<double> Tn = new List<double>();
            List<double> circTn = new List<double>();
            List<double> regenDelay = new List<double>();
            List<double> condPe = new List<double>();
            List<double> AccPe = new List<double>();
            List<double> Pn = new List<double>();
            List<double> index = new List<double>();
            List<double> Tnratio = new List<double>();

            List<double> Id = new List<double>();
            List<double> Iq = new List<double>();
            List<double> Pekwh = new List<double>();
            List<double> Pekws = new List<double>();
            List<double> Ploss = new List<double>();
            List<double> Pcu = new List<double>();
            List<double> Pfe = new List<double>();
            List<double> Pf = new List<double>();
            List<double> Pstr = new List<double>();
            List<double> Pw = new List<double>();
            List<double> Pinv = new List<double>();
            List<double> n = new List<double>();
            List<double> cdcI = new List<double>();
            List<double> Pbatt = new List<double>();
            List<double> AccPbatt = new List<double>();
            List<double> Temp = new List<double>();
            List<bool> nonNaNIndex = new List<bool>();

            //1
            double cumulativeSum = 0;
            double cumulativeSumPbtt = 0;

            double nacount = 0;
            double count = 0;
            double WLTCconpow=0;
            double WLTCdist=0;
            double bttcon = 0;
            double crusdist = 0;
            int lastNonNaNIndex;


            acckm.Add(ip["Vx"][0]);

            // calculate the differences and add them to acckm
            for (int i = 1; i < ip["Vx"].Count; i++)
            {
                acckm.Add(ip["Vx"][i] - ip["Vx"][i - 1]);
            }
            for (int i = 0; i < acckm.Count; i++)
            {
                acc.Add (acckm[i] / 3.6);
                dist.Add( ip["Vx"][i] * 1 / 3.6 * 1 + acc[i] * 1 / 2);
                wwrpm.Add(dist[i] / (2 * Math.PI * tire["Tirerw"]) * 60);
                ww.Add(wwrpm[i] / 60 * 2 * Math.PI);
                alpha.Add(acc[i] / tire["Tirerw"]);
                InertialL.Add(ip["GrossWt"] * acc[i]);
                Frolling.Add(0);
                if (ip["crsPtn"][i] > 0)
                {
                    Frolling[i] = ip["fr"] * ip["GrossWt"] * ip["gravity"] * slopecos;
                }
                //FAero = 0.5 * ip["p"] * ip["Af"] * ip["Cd"] * Math.Pow((ip["Vx"][i] + ip["Vwind"]) / 3.6, 2);
                FAero.Add( 0.5 * ip["p"] * ip["Af"] * ip["Cd"] * Math.Pow((ip["Vx"][i] + ip["Vwind"]) / 3.6, 2));
                Fx.Add(InertialL[i] + Fgravity + Frolling[i] + FAero[i]);

                if (ip["VDswitch"] == 1)
                {
                    Tw.Add(ip["J"] * alpha[i] + tire["Tirerw"] * (Fgravity + Frolling[i] + FAero[i]));
                }
                else 
                {
                    Tw.Add(ip["J"] * alpha[i]);
                }

                Pe.Add(ww[i] * Tw[i] / 1000 / 3600);
                wr.Add(ww[i]* gear["gdr"]);
                we.Add(wr[i] * pole / 2);
                werpm.Add(we[i] * 60 / ((pole / 2) * 2 * Math.PI));
                Tn.Add(Tw[i]/ gear["gdr"]);
                //Tn.CopyTo(0, circTn.ToArray(), 1, Tn.Count - 1);
                if (i > 0)
                {
                    circTn.Add( Tn[i - 1]);
                }
                else
                {
                    circTn.Add(0);
                }
                regenDelay.Add(0);
                if (ip["zone"] == 1)
                {
                    if (circTn[i] >= 0 && Tn[i] < 0)
                    {
                        regenDelay.Add(1);
                    }
                }
                //recheck condpe
                condPe.Add(-Pe[i]);
                if (ip["Vx"][i] > btt["regen_limit"])
                {
                    condPe.Add(-Pe[i] * (1 - btt["regen_ratio"]));
                }
                if (alpha[i] >= 0)
                {
                    condPe.Add(Pe[i]);
                }
                Pn.Add(1);
                if (Tn[i] < 0)
                {
                    Pn.Add(-1) ;
                }
                index.Add(0);
                cumulativeSum += -condPe[i];
                AccPe.Add(cumulativeSum);
                ipm.Tn.Add(0);
                for (int iter = 1; iter < 4100; iter++)
                {
                    //tmprpm[0] = 0;
                    //tmprpm.Add(ipm.rpm[iter] - ipm.increment);

                    if (rpmflag[iter - 1] == 1 &&
                        Math.Abs(werpm[i]) >= tmprpm[iter - 1] &&
                        Math.Abs(werpm[i]) < ipm.rpm[iter - 1] &&
                        Math.Abs(Tn[i]) >= ipm.Tn[iter - 1])
                    {
                        index[i] = iter;
                        //Console.WriteLine(index[i]);
                        // break;
                    }
                }
                ipm.rpm.RemoveAt(ipm.rpm.Count - 1);
                Tnratio.Add(0);
                Pekwh.Add(0);
                Pekws.Add(0);
                Id.Add(0);
                Iq.Add(0);
                Id.Add(0);
                Ploss.Add(0);
                Pcu.Add(0);
                Pfe.Add(0);
                Pstr.Add(0);
                Pf.Add(0);
                Pw.Add(0);
                Pinv.Add(0);
                n.Add(0);
                cdcI.Add(0);
                Pbatt.Add(0);
                Temp.Add(0);
                AccPbatt.Add(0);


                if (index[i] != 0)
                {
                    if (Tn[i] != 0)
                    {
                        Tnratio[i] = (Math.Abs(Tn[i]) - ipm.Tn[(int)index[i]-1]) / (ipm.Tn[(int)index[i]] - ipm.Tn[(int)index[i] - 1]);
                        //qDebug() << i << "WLTC error 1 : "<<index[i]; 
                        //qDebug() << i << "WLTC error 1 : "<<index[i]; 
                    }
                    Id[i] = (ipm.Id[(int)index[i]] - ipm.Id[(int)index[i] - 1]) * Tnratio[i] + ipm.Id[(int)index[i] - 1];
                    Iq[i] = Pn[i] * ((ipm.Iq[(int)index[i]] - ipm.Iq[(int)index[i] - 1]) * Tnratio[i] + ipm.Iq[(int)index[i - 1]]);
                    Pekwh[i] = Pn[i] * ((ipm.Pe[(int)index[i]] - ipm.Pe[(int)index[i] - 1]) * Tnratio[i] + ipm.Pe[(int)index[i - 1]]); // Col AL
                    Pekws[i] = Pekwh[i] / 3600; // Col AM
                    Ploss[i] = Pn[i] * ((ipm.Ploss[(int)index[i]] - ipm.Ploss[(int)index[i] - 1]) * Tnratio[i] + ipm.Ploss[(int)index[i - 1]]) / 3600; // Col AN
                    Pcu[i] = (ipm.Pcu[(int)index[i]] - ipm.Pcu[(int)index[i] - 1]) * Tnratio[i] + ipm.Pcu[(int)index[i - 1]]; // Col A0
                    Pfe[i] = (ipm.Pfe[(int)index[i]] - ipm.Pfe[(int)index[i] - 1]) * Tnratio[i] + ipm.Pfe[(int)index[i - 1]]; // Col AP
                    Pstr[i] = (ipm.Pstr[(int)index[i]] - ipm.Pstr[(int)index[i] - 1]) * Tnratio[i] + ipm.Pstr[(int)index[i - 1]]; // Col AQ
                    Pf[i] = (ipm.Pfric[(int)index[i]] - ipm.Pfric[(int)index[i] - 1]) * Tnratio[i] + ipm.Pfric[(int)index[i - 1]]; // Col AR
                    Pw[i] = (ipm.Pwind[(int)index[i]] - ipm.Pwind[(int)index[i] - 1]) * Tnratio[i] + ipm.Pwind[(int)index[i - 1]]; // Col AS
                    Pinv[i] = (ipm.Pinv[(int)index[i]] - ipm.Pinv[(int)index[i] - 1]) * Tnratio[i] + ipm.Pinv[(int)index[i - 1]]; // Col AT

                    n[i] = (ipm.posin[(int)index[i]] - ipm.posin[(int)index[i] - 1]) * Tnratio[i] + ipm.posin[(int)index[i - 1]]; // Col AX
                    cdcI[i] = Math.Sqrt(Math.Pow((Id[i]), 2) + Math.Pow((Iq[i]), 2)); // Col AV

                    if (Pn[i] > 0)
                    {
                        Pbatt[i] = Pekws[i] / n[i];
                    }
                    else if(ip["Vx"][i] > btt["regen_limit"] && regenDelay[i] != 1)
                    {
                        Pbatt[i] = -(Pekws[i] / n[i] - (Pekws[i] - Math.Abs(Ploss[i])) * btt["regen_ratio"]);
                    }
                    else
                    {
                        Pbatt[i] = -Pekws[i] / n[i];
                    }
                    
                }
                else
                {
                    if (Tn[i] < 0)
                    {
                        nacount += 1; 
                    }
                }
                if(i>0) {
                    if (n[i] > 0)
                    {
                        Temp[i] = (ipmclass.Temp[(int)index[i] + 1] - ipmclass.Temp[(int)index[i]]) * Tnratio[i] + ipmclass.Temp[(int)index[i]];
                        count = count + 1;
                    }
                    else
                    {
                        Temp[i] = Temp[i - 1];
                    }
                    AccPbatt[i] = AccPbatt[i - 1] - Pbatt[i];
                }
                else
                {
                    AccPbatt[i] = 0;
                    Temp[i] = 57.92;
                }
                
               // cumulativeSumPbtt -= Pbatt[i];
               // AccPbatt[i] = cumulativeSumPbtt;
            //nonNaNIndex = AccPbatt.Select(value => !double.IsNaN(value)).ToList();
            //lastNonNaNIndex = nonNaNIndex.LastIndexOf(true);
            //if (lastNonNaNIndex != -1) // Check if there's at least one non-NaN value
            //{
            //    WLTCconpow = -AccPbatt[lastNonNaNIndex];
            //}
            }
            //WLTCdist = dist.Where(d => !double.IsNaN(d)).Sum() / 1000;
            //bttcon = WLTCconpow / WLTCdist * 10;
            //crusdist = btt["charge"] / bttcon * 10;
            //Console.WriteLine(bttcon);
            WLTCconpow = -AccPbatt.Last(); // Using Last() to get the last element from the list and negate it
            WLTCdist = dist.Sum() / 1000; // Sum up the dist list and divide by 1000
            bttcon = WLTCconpow / WLTCdist * 10; // Calculate bttcon
            crusdist = btt["charge"] / bttcon * 10; //
            Console.WriteLine(bttcon);





        }
        public void constTrq(List<double[]> arr, ipmclass ipm)
        {
            List<int> ctwerpm = new List<int>(Enumerable.Range(1000, 19000).Where(x => x % 1000 == 0))
                        .Concat(Enumerable.Range(1000, 19000).Where(x => x % 1000 == 0))
                        .Concat(Enumerable.Range(500, 6501).Where(x => x % 500 == 0))
                        .Concat(Enumerable.Range(7100, 401).Where(x => x % 100 == 0))
                        .Concat(Enumerable.Range(500, 4001).Where(x => x % 500 == 0))
                        .Concat(Enumerable.Range(4520, 181).Where(x => x % 20 == 0))
                        .ToList();

            List<double> ctTn = new List<double>(new double[76]);
            ctTn.InsertRange(0, arr[0].Take(19));    // Copy arr[0] to ctTn[0:19]
            ctTn.InsertRange(19, arr[1].Take(19));   // Copy arr[1] to ctTn[19:38]
            ctTn.InsertRange(38, arr[2].Take(19));   // Copy arr[2] to ctTn[38:57]
            ctTn.InsertRange(57, arr[3].Take(19));   // Copy arr[3] to ctTn[57:76]

            List<int> ctindex = new List<int>(new int[ctwerpm.Count]);
            List<double> cTratio = new List<double>(new double[ctwerpm.Count]);

            // Initialize ctId and ctIq as lists
            List<double> ctId = new List<double>(new double[ctwerpm.Count]);
            List<double> ctIq = new List<double>(new double[ctwerpm.Count]);

            // Initialize tmprpm and rpmflag lists with a size of 4100
            List<double> tmprpm = new List<double>(new double[4100]);
            List<int> rpmflag = new List<int>(Enumerable.Repeat(1, 4100));
            for (int idx = 0; idx < ctwerpm.Count; idx++)
            {
                // Append an extra element to ipm.Tn
                ipm.Tn.Add(0);  // Assuming ipm.Tn is a List<double>

                // Iterate through RPMs
                for (int iter = 1; iter < 4100; iter++)
                {
                    tmprpm[0] = 0;
                    tmprpm[iter] = ipm.rpm[iter] - ipm.increment;

                    if (rpmflag[iter - 1] == 1 &&
                        Math.Abs(ctwerpm[idx]) >= tmprpm[iter - 1] &&
                        Math.Abs(ctwerpm[idx]) < ipm.rpm[iter - 1] &&
                        Math.Abs(ctTn[idx]) >= ipm.Tn[iter - 1])
                    {
                        ctindex[idx] = iter;
                    }
                }

                // Remove the extra element from ipm.Tn
                ipm.Tn.RemoveAt(ipm.Tn.Count - 1);

                // Perform calculations
                if (ctindex[idx] != 0)
                {
                    if (ctTn[idx] != 0)
                    {
                        cTratio[idx] = (Math.Abs(ctTn[idx]) - ipm.Tn[ctindex[idx] - 1]) /
                                       (ipm.Tn[ctindex[idx]] - ipm.Tn[ctindex[idx] - 1]);

                        ctId[idx] = (ipm.Id[ctindex[idx]] - ipm.Id[ctindex[idx] - 1]) * cTratio[idx] + ipm.Id[ctindex[idx] - 1];
                        ctIq[idx] = (ipm.Iq[ctindex[idx]] - ipm.Iq[ctindex[idx] - 1]) * cTratio[idx] + ipm.Iq[ctindex[idx] - 1];
                    }
                }
                else
                {
                    ctId[idx] = double.NaN;  // Use NaN to represent None
                    ctIq[idx] = double.NaN;
                }
            }
        }






    }
}
