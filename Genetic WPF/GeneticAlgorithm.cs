using ExcelDataReader;
using System;
using System.Collections.Generic;
using System.Data;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection.PortableExecutable;
using System.Security.Cryptography;
using System.Security.Policy;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Windows.Controls;
using System.Windows.Shapes;
using System.Windows.Xps.Packaging;
using static System.Runtime.InteropServices.JavaScript.JSType;

namespace Genetic_WPF
{

    class GeneticAlgorithm
    {
        private MainWindow _mainWindow;
        public static Dictionary<string, double> mot = new Dictionary<string, double>();
        public static Dictionary<string, object> dynamic_variables = new Dictionary<string, object>();
        public static Dictionary<string, object> dynamic_variables1 = new Dictionary<string, object>();
        public static List<string> unticked_titles = new List<string>();
        public static List<string> ticked_titles = new List<string>();
        public static Dictionary<string, (double Min, double Max)> slider_bounds = new Dictionary<string, (double, double)>();
        public static Dictionary<string, (double Min, double Max)> slider_bounds1 = new Dictionary<string, (double, double)>();
        public static Dictionary<string, (double Min, double Max)> slider_bounds2 = new Dictionary<string, (double, double)>();
        public static Dictionary<string, bool> ticked_values = new Dictionary<string, bool>();
        public static Dictionary<string, object> individual_values = new Dictionary<string, object>();
        public static Dictionary<string, object> unticked_values = new Dictionary<string, object>();
        public static Dictionary<string, object> other_values = new Dictionary<string, object>();
        public static Dictionary<string, bool> entry_checkbox_map = new Dictionary<string, bool>();
        public static Dictionary<string, (double Min, double Max)> slider_bounds_obj = new Dictionary<string, (double, double)>();
        public static Dictionary<string, (double Min, double Max)> slider_bounds_wltc_obj = new Dictionary<string, (double, double)>();
        public static Dictionary<string, (double Min, double Max)> slider_bounds_wltc = new Dictionary<string, (double, double)>();
        public static int clicks = 0;
        //public static MotorWindow motorWindow;
        public static double selected_slider_value = 0;
        public static Dictionary<string, object> dynamic_variables_wltc = new Dictionary<string, object>();
        public static Dictionary<string, object> dynamic_variables_wltc1 = new Dictionary<string, object>();
        public static Dictionary<string, object> Output_result_wltc = new Dictionary<string, object>();
        public static List<string> unticked_titles_wltc = new List<string>();
        public static List<string> ticked_titles_wltc = new List<string>();
        public static Dictionary<string, bool> ticked_values_wltc = new Dictionary<string, bool>();
        public static Dictionary<string, object> individual_values_wltc = new Dictionary<string, object>();
        public static Dictionary<string, object> unticked_values_wltc = new Dictionary<string, object>();
        public static Dictionary<string, object> other_values_wltc = new Dictionary<string, object>();
        public static Dictionary<string, bool> entry_checkbox_map_wltc = new Dictionary<string, bool>();
        public static int clicks_wltc = 0;

        //public static Dictionary<string, double>? mot;
        //public static Dictionary<string, double>? dynamic_variables;
        //public static Dictionary<string, double>? dynamic_variables1;
        //UI components with static values
        public double? EffMaxSpeed { get; set; }  //double

        public GeneticAlgorithm(MainWindow mainWindow)
        {
            _mainWindow = mainWindow;

            //Dictionary<string, double> mot = new Dictionary<string, double>();
            //Dictionary<string, double> mot = new Dictionary<string, double>();
            //Dictionary<string, object> dynamic_variables = new Dictionary<string, object>();
            //Dictionary<string, object> dynamic_variables1 = new Dictionary<string, object>();
            //List<string> unticked_titles = new List<string>();
            //List<string> ticked_titles = new List<string>();
            //Dictionary<string, (double Min, double Max)> slider_bounds = new Dictionary<string, (double, double)>();
            //Dictionary<string, (double Min, double Max)> slider_bounds1 = new Dictionary<string, (double, double)>();
            //Dictionary<string, (double Min, double Max)> slider_bounds2 = new Dictionary<string, (double, double)>();
            //Dictionary<string, bool> ticked_values = new Dictionary<string, bool>();
            //Dictionary<string, object> individual_values = new Dictionary<string, object>();
            //Dictionary<string, object> unticked_values = new Dictionary<string, object>();
            //Dictionary<string, object> other_values = new Dictionary<string, object>();
            //Dictionary<string, bool> entry_checkbox_map = new Dictionary<string, bool>();
            //Dictionary<string, (double Min, double Max)> slider_bounds_obj = new Dictionary<string, (double, double)>();
            //Dictionary<string, (double Min, double Max)> slider_bounds_wltc_obj = new Dictionary<string, (double, double)>();
            //Dictionary<string, (double Min, double Max)> slider_bounds_wltc = new Dictionary<string, (double, double)>();
            //int clicks = 0;
            ////public static MotorWindow motorWindow;
            //double selected_slider_value = 0;
            //Dictionary<string, object> dynamic_variables_wltc = new Dictionary<string, object>();
            //Dictionary<string, object> dynamic_variables_wltc1 = new Dictionary<string, object>();
            //Dictionary<string, object> Output_result_wltc = new Dictionary<string, object>();
            //List<string> unticked_titles_wltc = new List<string>();
            //List<string> ticked_titles_wltc = new List<string>();
            //Dictionary<string, bool> ticked_values_wltc = new Dictionary<string, bool>();
            //Dictionary<string, object> individual_values_wltc = new Dictionary<string, object>();
            //Dictionary<string, object> unticked_values_wltc = new Dictionary<string, object>();
            //Dictionary<string, object> other_values_wltc = new Dictionary<string, object>();
            //Dictionary<string, bool> entry_checkbox_map_wltc = new Dictionary<string, bool>();
            //int clicks_wltc = 0;
        }

        public static DataTable ReadExcelSheet(string filePath, string sheetName)
        {
            // Open the Excel file
            using (var stream = File.Open(filePath, FileMode.Open, FileAccess.Read))
            {
                // Auto-detect format, supports .xls, .xlsx, etc.
                using (var reader = ExcelReaderFactory.CreateReader(stream))
                {
                    // Configure to read the Excel file into a DataSet (multiple sheets supported)
                    var result = reader.AsDataSet(new ExcelDataSetConfiguration()
                    {
                        ConfigureDataTable = (_) => new ExcelDataTableConfiguration()
                        {
                            UseHeaderRow = true // Use the first row as column headers
                        }
                    });

                    // Get the specified sheet by name
                    DataTable dataTable = result.Tables[sheetName];
                    return dataTable;
                }
            }
        }

        private string DataTableToString(DataTable dataTable)
        {
            StringBuilder sb = new StringBuilder();

            // Add column headers
            foreach (DataColumn column in dataTable.Columns)
            {
                sb.Append(column.ColumnName + "\t");
            }
            sb.AppendLine();

            // Add rows
            foreach (DataRow row in dataTable.Rows)
            {
                foreach (var item in row.ItemArray)
                {
                    sb.Append(item + "\t");
                }
                sb.AppendLine();
            }

            return sb.ToString();
        }


        public static string Convert1DListToString(List<double> list1D)
        {
            StringBuilder sb = new StringBuilder();

            int rows = list1D.Count; // Number of rows

            // Iterate over the 2D array and append the values to StringBuilder
            for (int i = 0; i < rows; i++)
            {
                sb.Append(list1D[i]); // Use tab to separate columns
                sb.AppendLine(); // Move to the next line after each row
            }

            return sb.ToString();
        }

        public string Convert2DListToString(List<List<double>> list2D)
        {
            StringBuilder sb = new StringBuilder();

            int rows = list2D.Count; // Number of rows
            int columns = list2D[0].Count; // Number of columns


            // Iterate over the 2D array and append the values to StringBuilder
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    sb.Append(list2D[i][j] + "\t"); // Use tab to separate columns
                }
                sb.AppendLine(); // Move to the next line after each row
            }

            return sb.ToString();
        }

        public string Convert3DListToString(List<List<List<double>>> list3D)
        {
            StringBuilder sb = new StringBuilder();

            int rows = list3D.Count; // Number of rows
            int columns = list3D[0].Count; // Number of columns
            int columns1 = list3D[0][0].Count; // Number of columns


            // Iterate over the 2D array and append the values to StringBuilder
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    for (int K = 0; K < columns1; K++)
                        sb.Append(list3D[i][j][K] + "\t"); // Use tab to separate columns
                }
                sb.AppendLine(); // Move to the next line after each row
            }

            return sb.ToString();
        }


        private string Convert2DArrayToString(double[,] array2D)
        {
            StringBuilder sb = new StringBuilder();

            int rows = array2D.GetLength(0); // Number of rows
            int columns = array2D.GetLength(1); // Number of columns

            // Iterate over the 2D array and append the values to StringBuilder
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    sb.Append(array2D[i, j] + "\t"); // Use tab to separate columns
                }
                sb.AppendLine(); // Move to the next line after each row
            }

            return sb.ToString();
        }

        private string Convert1DArrayToString(object[] array)
        {
            StringBuilder sb = new StringBuilder();

            int rows = array.GetLength(0); // Number of rows

            for (int i = 0; i < rows; i++)
            {
                sb.Append(array[i]); // Append each element followed by a newline
                sb.AppendLine(); // Append each element followed by a newline

            }

            return sb.ToString();
        }

        private string ConvertDictionaryToString(Dictionary<string, object> dictionary)
        {
            StringBuilder sb = new StringBuilder();

            foreach (var kvp in dictionary)
            {
                sb.AppendLine($"{kvp.Key}:");

                // Assuming each value is an array, convert the array to a string
                if (kvp.Value is Array array)
                {
                    foreach (var item in array)
                    {
                        sb.AppendLine($"\t{item}");
                    }
                }
                else
                {
                    sb.AppendLine($"\t{kvp.Value}");
                }

                sb.AppendLine(); // New line after each dictionary entry
            }

            return sb.ToString();
        }

        private string ConvertDictionaryKeyToString(Dictionary<string, object> dictionary, string key)
        {
            // Check if the dictionary contains the key
            if (dictionary.ContainsKey(key))
            {
                // Retrieve the value associated with the key
                var value = dictionary[key];

                // Check if the value is an array
                if (value is Array array)
                {
                    // Convert the array to a string with each element on a new line
                    StringBuilder sb = new StringBuilder();
                    foreach (var item in array)
                    {
                        sb.AppendLine(item?.ToString());
                        //sb.AppendLine(item);

                    }
                    return sb.ToString();
                }
                else
                {
                    // If it's not an array, just return the value as a string
                    return value.ToString();
                }
            }
            else
            {
                return $"Key '{key}' not found in the dictionary.";
            }
        }

        public List<List<double>> DigiTwin(double pole, double Vm, List<List<double>> DTtab) {
            //DTtab = np.array(DTtab)
            double[] brr = { 0.968309134210145, 1.10411600242257, 0.000289328513383774, 0.000382804854148081, 0.126188070323104, 0.000000 };
            double[] arr = brr;

            arr[0] = arr[0] / ((3.0 / 2.0) * pole / 2);
            arr[1] = arr[1] / ((3.0 / 2.0) * pole / 2);
            arr[2] = arr[2] * 1000000;
            arr[3] = arr[3] * 1000000;

            double Ld2 = 0.000212744959421379;
            double Lq2 = 0.000311235559094195;
            double psi2 = 0.0729606243955028;

            // Update digital twin table
            for (int i = 0; i < 6; i++)
                DTtab[2][i] = arr[i];
            // print(DTtab)
            DTtab[5][2] = Ld2 * 1000000;
            DTtab[5][3] = Lq2 * 1000000;
            DTtab[5][4] = psi2;

            return DTtab;
        }


        public Dictionary<string, List<List<List<double>>>> PsiLdLq(double pole, List<List<double>> DTtab, int Hdcalf, int Hqcalf, int type, int DTnum, double AdjPM, double AdjLd, double AdjLq, List<List<double>> Psid, List<List<double>> Psiq, List<List<double>> Omega)

        {
            Dictionary<string, List<List<List<double>>>> PLmat = new Dictionary<string, List<List<List<double>>>>();
            
            List<double> _Idarr = new List<double>();
            for (int i = -650; i < 1; i+=25)
                _Idarr.Add(i);
            double[] Idarr = _Idarr.ToArray();

            List<double> _Iqarr = new List<double>();
            for (int i = 650; i > -1; i-=25)
                _Iqarr.Add(i);

            double[] Iqarr = _Iqarr.ToArray();
            double ndrr = Idarr.GetLength(0);
            double nqrr = Iqarr.GetLength(0);

            AdjPM = AdjPM / 100.0;
            AdjLd = AdjLd / 100.0;
            AdjLq = AdjLq / 100.0;

            List<List<double>> Hdmat = new List<List<double>>(Enumerable.Range(0, 6).Select(x => new List<double>(new double[Idarr.GetLength(0)])));

            for (int ind = 0; ind < Idarr.GetLength(0); ind++) 
            {
                for (int idx = 0; idx < 6; idx++)
                    Hdmat[idx][ind] = DTtab[idx][4] + 0.000001 * DTtab[idx][2] * Idarr[ind]; 
                if (Hdcalf == 2)
                    Hdmat[5][ind] = 5.375 * Math.Pow(10, -8) * Math.Pow(Idarr[ind], 2) + 2.11 * Math.Pow(10, -4) * Idarr[ind] + 7.305 * Math.Pow(10, -2);
            }

            List<List<double>> Hqmat = new List<List<double>>(Enumerable.Range(0, 6).Select(x => new List<double>(new double[Idarr.GetLength(0)])));

            for (int ind = 0; ind < Iqarr.GetLength(0); ind++)
            {
                for (int idx = 0; idx < 6; idx++)
                    Hqmat[idx][ind] = 0.000001 * DTtab[idx][3] * Iqarr[ind] + DTtab[idx][5];
                if (Hqcalf == 2)
                    Hqmat[5][ind] = 4.903 * Math.Pow(10, -13) * Math.Pow(Iqarr[ind], 3) + 1.468 * Math.Pow(10, -6) * Math.Pow(Iqarr[ind], 2) + 6.678 * Math.Pow(10, -4) * Iqarr[ind];
            }

            double solverA = 0.845762136019501;
            double solverB = 1.37338000111102;

            //    # Assuming numel(Iqarr) and numel(Idarr) are defined before this code block
            List<List<double>> H2onemat = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> H2twomat = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> H2threemat = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> H2fourmat = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> H2fivemat = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> H2sixmat = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> H2sevenmat = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));

            List<List<double>> oneK0 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> twoK0 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> threeK0 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> fourK0 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> fiveK0 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> sixK0 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> sevenK0 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));


            //_mainWindow.myTextbox.Text = Psiq.Count.ToString();
            //_mainWindow.myTextbox.Text += "," + Psiq[0].Count.ToString();

            //_mainWindow.myTextbox.Text = Iqarr.GetLength(0).ToString();
            //_mainWindow.myTextbox.Text += "," + Idarr.GetLength(0).ToString();

            for (int ind = 0; ind < Idarr.GetLength(0); ind++)
                for (int idx = 0; idx < Iqarr.GetLength(0); idx++)
                {
                    H2onemat[idx][ind] = Math.Pow((Hdmat[0][ind] / DTtab[0][0]), 2) + Math.Pow((Hqmat[0][idx] / DTtab[0][1]), 2);
                    H2twomat[idx][ind] = Math.Pow((Hdmat[1][ind] / DTtab[1][0]), 2) + Math.Pow((Hqmat[1][idx] / DTtab[1][1]), 2);
                    H2threemat[idx][ind] = Math.Pow((Hdmat[2][ind] / DTtab[2][0]), 2) + Math.Pow((Hqmat[2][idx] / DTtab[2][1]), 2);
                    H2fourmat[idx][ind] = Math.Pow((Hdmat[3][ind] / DTtab[3][0]), 2) + Math.Pow((Hqmat[3][idx] / DTtab[3][1]), 2);
                    H2fivemat[idx][ind] = Math.Pow((Hdmat[4][ind] / DTtab[4][0]), 2) + Math.Pow((Hqmat[4][idx] / DTtab[4][1]), 2);
                    H2sixmat[idx][ind] = Omega[idx][ind] * Math.Pow((pole / 2.0 * 3.0 / 2.0), 2);
                    H2sevenmat[idx][ind] = Math.Pow((pole / 2.0 * 3.0 / 2.0), 2.0) * (Math.Pow((Hdmat[5][ind] / solverA), 2) + Math.Pow((Hqmat[5][idx] / solverB), 2));

                    oneK0[idx][ind] = Math.Sqrt(H2onemat[idx][ind]);
                    twoK0[idx][ind] = Math.Sqrt(H2twomat[idx][ind]);
                    threeK0[idx][ind] = Math.Sqrt(H2threemat[idx][ind]);
                    fourK0[idx][ind] = Math.Sqrt(H2fourmat[idx][ind]);
                    fiveK0[idx][ind] = Math.Sqrt(H2fivemat[idx][ind]);
                    sixK0[idx][ind] = Math.Sqrt(H2sixmat[idx][ind]);
                    sevenK0[idx][ind] = Math.Sqrt(H2sevenmat[idx][ind]);
                }


            DataTable _aval = ReadExcelSheet("InputTableFile.xlsx", "a(k)");

            ////var Psidstring = DataTableToString(_Psid);
            ////_mainWindow.myTextbox.Text = Psidstring;

            List<List<double>> aval = new List<List<double>>(Enumerable.Range(0, _aval.Rows.Count).Select(x => new List<double>(new double[_aval.Columns.Count])));

            for (int i = 0; i < _aval.Rows.Count; i++)  //Skip one row
                for (int j = 0; j < _aval.Columns.Count; j++)   //Skip one column
                {
                    aval[i][j] = (double)_aval.Rows[i][j];
                }

            //var testString = Convert2DArrayToString(aval);
            //_mainWindow.myTextbox.Text = testString;

            List<List<double>> oneKmat = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> Kk0one = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> twoKmat = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> Kk0two = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> threeKmat = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> Kk0three = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> fourKmat = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> Kk0four = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> fiveKmat = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> Kk0five = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> sixKmat = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> Kk0six = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> sevenKmat = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> Kk0seven = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> eightKmat = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> Kk0eight = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));

            for (int idx = 0; idx < Iqarr.GetLength(0); ++idx)
            {
                for (int ind = 0; ind < Idarr.GetLength(0); ++ind)
                {
                    oneKmat[idx][ind] = DTtab[0][6] * oneK0[idx][ind] / Math.Sqrt(1 + Math.Pow(DTtab[0][6] * oneK0[idx][ind], 2));
                    Kk0one[idx][ind] = oneKmat[idx][ind] / oneK0[idx][ind];

                    twoKmat[idx][ind] = Math.Tanh(DTtab[1][6] * twoK0[idx][ind]);
                    Kk0two[idx][ind] = twoKmat[idx][ind] / twoK0[idx][ind];

                    threeKmat[idx][ind] = DTtab[2][6] * threeK0[idx][ind] / Math.Sqrt(1 + Math.Pow(DTtab[2][6] * threeK0[idx][ind], 2));
                    Kk0three[idx][ind] = threeKmat[idx][ind] / threeK0[idx][ind];

                    fourKmat[idx][ind] = Math.Atan(DTtab[3][6] * fourK0[idx][ind]);
                    Kk0four[idx][ind] = fourKmat[idx][ind] / fourK0[idx][ind];

                    fiveKmat[idx][ind] = 1 - Math.Exp(-DTtab[4][6] * fiveK0[idx][ind]);
                    Kk0five[idx][ind] = fiveKmat[idx][ind] / fiveK0[idx][ind];

                    sixKmat[idx][ind] = (1.32962896540908) * sixK0[idx][ind] / Math.Sqrt(1 + Math.Pow((1.32962896540908) * sixK0[idx][ind], 2));
                    Kk0six[idx][ind] = sixKmat[idx][ind] / sixK0[idx][ind];

                    sevenKmat[idx][ind] = DTtab[5][6] * sevenK0[idx][ind] / Math.Sqrt(1 + Math.Pow(DTtab[5][6] * sevenK0[idx][ind], 2));
                    Kk0seven[idx][ind] = sevenKmat[idx][ind] / sevenK0[idx][ind];

                    eightKmat[idx][ind] = Math.Atan(aval[idx][0] * sevenK0[idx][ind]); // Assuming aval is properly defined
                    Kk0eight[idx][ind] = eightKmat[idx][ind] / sevenK0[idx][ind];
                }
            }

            //var testString = Convert2DListToString(eightKmat);
            //_mainWindow.myTextbox.Text = testString;

            List<List<double>> Tlin1 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> Tlin2 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> Tlin3 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> Tlin4 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> Tlin5 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> Tlin7 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> TKk0cal6 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));

            for (int ind = 0; ind < Idarr.GetLength(0); ++ind)
            {
                for (int idx = 0; idx < Iqarr.GetLength(0); ++idx)
                {
                    Tlin1[idx][ind] = (3.0 / 2.0) * (pole / 2.0) * (Hdmat[0][ind] * Iqarr[idx] - Hqmat[0][idx] * Idarr[ind]);
                    Tlin2[idx][ind] = (3.0 / 2.0) * (pole / 2.0) * (Hdmat[1][ind] * Iqarr[idx] - Hqmat[1][idx] * Idarr[ind]);
                    Tlin3[idx][ind] = (3.0 / 2.0) * (pole / 2.0) * (Hdmat[2][ind] * Iqarr[idx] - Hqmat[2][idx] * Idarr[ind]);
                    Tlin4[idx][ind] = (3.0 / 2.0) * (pole / 2.0) * (Hdmat[3][ind] * Iqarr[idx] - Hqmat[3][idx] * Idarr[ind]);
                    Tlin5[idx][ind] = (3.0 / 2.0) * (pole / 2.0) * (Hdmat[4][ind] * Iqarr[idx] - Hqmat[4][idx] * Idarr[ind]);
                    Tlin7[idx][ind] = (3.0 / 2.0) * (pole / 2.0) * (Hdmat[5][ind] * Iqarr[idx] - Hqmat[5][idx] * Idarr[ind]);
                    TKk0cal6[ind][idx] = (3.0 / 2.0) * (pole / 2.0) * (Iqarr[idx] * Psid[ind][idx] - Idarr[idx] * Psiq[ind][idx]);
                }
            }

            DataTable _t_ave_df = ReadExcelSheet("InputTableFile.xlsx", "T_ave_matrix");

            ////var Psidstring = DataTableToString(_Psid);
            ////_mainWindow.myTextbox.Text = Psidstring;

            List<List<double>> t_ave = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));

            for (int i = 1; i < _t_ave_df.Rows.Count; i++)  //Skip one row
                for (int j = 2; j < _t_ave_df.Columns.Count; j++)   //Skip one column
                {
                    t_ave[i-1][j-2] = (double)_t_ave_df.Rows[i][j];
                }

            //string testString = t_ave.Count.ToString()+"x"+ t_ave[0].Count.ToString();
            //_mainWindow.myTextbox.Text = testString;
            //testString = Iqarr.GetLength(0).ToString() + "x" + Idarr.GetLength(0).ToString();
            //_mainWindow.myTextbox.Text += testString;

            List<List<double>> TKk0cal1 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> TKk0cal2 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> TKk0cal3 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> TKk0cal4 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> TKk0cal5 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> TKk0cal7 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> TKk0cal8 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));

            for (int idx = 0; idx < Iqarr.GetLength(0); ++idx)
            {
                for (int ind = 0; ind < Idarr.GetLength(0); ++ind)
                {
                    TKk0cal1[idx][ind] = Tlin1[idx][ind] * Kk0one[idx][ind];
                    TKk0cal2[idx][ind] = Tlin2[idx][ind] * Kk0two[idx][ind];
                    TKk0cal3[idx][ind] = Tlin3[idx][ind] * Kk0three[idx][ind];
                    TKk0cal4[idx][ind] = Tlin4[idx][ind] * Kk0four[idx][ind];
                    TKk0cal5[idx][ind] = Tlin5[idx][ind] * Kk0five[idx][ind];
                    TKk0cal7[idx][ind] = Tlin7[idx][ind] * Kk0seven[idx][ind];
                    TKk0cal8[idx][ind] = Tlin7[idx][ind] * Kk0eight[idx][ind]; // Note: Using Tlin7 instead of Tlin8 as per your Python code
                }
            }

            List<List<List<double>>> DTKk0 = new List<List<List<double>>>(Enumerable.Range(0, 8).Select(x => new List<List<double>>(Enumerable.Range(0, t_ave.Count).Select(y => new List<double>(new double[t_ave[0].Count])))));
            //std::vector<std::vector<std::vector<double>>> DTKk0(8, std::vector<std::vector<double>>(t_ave.size(), std::vector<double>(t_ave[0].size(), 0.0)));

            for (int j = 0; j < t_ave.Count; ++j)
            {
                for (int k = 0; k < t_ave[j].Count; ++k)
                {
                    DTKk0[0][j][k] = TKk0cal1[j][k] - t_ave[j][k];
                    DTKk0[1][j][k] = TKk0cal2[j][k] - t_ave[j][k];
                    DTKk0[2][j][k] = TKk0cal3[j][k] - t_ave[j][k];
                    DTKk0[3][j][k] = TKk0cal4[j][k] - t_ave[j][k];
                    DTKk0[4][j][k] = TKk0cal5[j][k] - t_ave[j][k];
                    DTKk0[5][j][k] = TKk0cal6[j][k] - t_ave[j][k];
                    DTKk0[6][j][k] = TKk0cal7[j][k] - t_ave[j][k];
                    DTKk0[7][j][k] = TKk0cal8[j][k] - t_ave[j][k];
                    // You can similarly assign values for other indices i = 1 to 7
                }
            }

            List<List<double>> xPsimOne = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> xPsimTwo = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> xPsimThree = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> xPsimFour = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> xPsimFive = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> xPsimSix = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> xPsimSeven = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));

            List<List<double>> xLd1 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> xLd2 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> xLd3 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> xLd4 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> xLd5 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> xLd6 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> xLd7 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));

            List<List<double>> xLq1 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> xLq2 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> xLq3 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> xLq4 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> xLq5 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> xLq6 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));
            List<List<double>> xLq7 = new List<List<double>>(Enumerable.Range(0, Iqarr.GetLength(0)).Select(x => new List<double>(new double[Idarr.GetLength(0)])));


            for (int ind = 0; ind < Idarr.GetLength(0); ++ind)
            {
                for (int idx = 0; idx < Iqarr.GetLength(0); ++idx)
                {
                    // Calculate xPsim values
                    xPsimOne[ind][idx] = Kk0one[ind][idx] * DTtab[0][4];
                    xPsimTwo[ind][idx] = Kk0two[ind][idx] * DTtab[1][4];
                    xPsimThree[ind][idx] = Kk0three[ind][idx] * DTtab[2][4];
                    xPsimFour[ind][idx] = Kk0four[ind][idx] * DTtab[3][4];
                    xPsimFive[ind][idx] = Kk0five[ind][idx] * DTtab[4][4];
                    xPsimSix[ind][idx] = Kk0six[ind][idx] * DTtab[4][4];
                    xPsimSeven[ind][idx] = Kk0seven[ind][idx] * DTtab[4][4];

                    // Calculate xLd values
                    xLd1[ind][idx] = Kk0one[ind][idx] * DTtab[0][2] * 0.000001;
                    xLd2[ind][idx] = Kk0two[ind][idx] * DTtab[1][2] * 0.000001;
                    xLd3[ind][idx] = Kk0three[ind][idx] * DTtab[2][2] * 0.000001;
                    xLd4[ind][idx] = Kk0four[ind][idx] * DTtab[3][2] * 0.000001;
                    xLd5[ind][idx] = Kk0five[ind][idx] * DTtab[0][2] * 0.000001;
                    xLd6[ind][idx] = Kk0six[ind][idx] * DTtab[5][2] * 0.000001;
                    xLd7[ind][idx] = Kk0seven[ind][idx] * DTtab[5][2] * 0.000001;

                    // Calculate xLq values
                    xLq1[ind][idx] = Kk0one[ind][idx] * DTtab[0][3] * 0.000001;
                    xLq2[ind][idx] = Kk0two[ind][idx] * DTtab[1][3] * 0.000001;
                    xLq3[ind][idx] = Kk0three[ind][idx] * DTtab[2][3] * 0.000001;
                    xLq4[ind][idx] = Kk0four[ind][idx] * DTtab[3][3] * 0.000001;
                    xLq5[ind][idx] = Kk0five[ind][idx] * DTtab[4][3] * 0.000001;
                    xLq6[ind][idx] = Kk0six[ind][idx] * DTtab[5][3] * 0.000001;
                    xLq7[ind][idx] = Kk0seven[ind][idx] * DTtab[5][3] * 0.000001;
                }
            }

            List<List<double>> xPsimOnearr = new List<List<double>>(Enumerable.Range(0, xPsimOne.Count).Select(x => new List<double>(new double[xPsimOne[0].Count])));
            List<List<double>> xLd1arr = new List<List<double>>(Enumerable.Range(0, xLd1.Count).Select(x => new List<double>(new double[xLd1[0].Count])));
            List<List<double>> xLq1arr = new List<List<double>>(Enumerable.Range(0, xLq1.Count).Select(x => new List<double>(new double[xLq1[0].Count])));


            for (int i = 0; i < xPsimOne.Count; ++i)
            {
                for (int j = 0; j < xPsimOne[0].Count; ++j)
                {
                    xPsimOnearr[i][j] = xPsimOne[i][j];
                    xLd1arr[i][j] = xLd1[i][j];
                    xLq1arr[i][j] = xLq1[i][j];
                }
            }


            List<List<double>> xPsiMode = new List<List<double>>(Enumerable.Range(0, 7).Select(x => new List<double>(new double[xPsimOnearr[0].Count])));
            List<List<double>> xLdMode = new List<List<double>>(Enumerable.Range(0, 7).Select(x => new List<double>(new double[xLd1arr[0].Count])));
            List<List<double>> xLqMode = new List<List<double>>(Enumerable.Range(0, 7).Select(x => new List<double>(new double[xLq1arr[0].Count])));

            if (type == 1)//xPsimOne[i].Average()
            {
                for (int i = 0; i < xPsimOne.Count; ++i)
                {
                    xPsiMode[0][i] = xPsimOne[i].Average();
                    xPsiMode[1][i] = xPsimTwo[i].Average();
                    xPsiMode[2][i] = xPsimThree[i].Average();
                    xPsiMode[3][i] = xPsimFour[i].Average();
                    xPsiMode[4][i] = xPsimFive[i].Average();
                    xPsiMode[5][i] = xPsimSix[i].Average();
                    xPsiMode[6][i] = xPsimSeven[i].Average();
                }

                // Calculate means for xLdMode
                for (int i = 0; i < xLd1.Count; ++i)
                {
                    xLdMode[0][i] = xLd1[i].Average();
                    xLdMode[1][i] = xLd2[i].Average();
                    xLdMode[2][i] = xLd3[i].Average();
                    xLdMode[3][i] = xLd4[i].Average();
                    xLdMode[4][i] = xLd5[i].Average();
                    xLdMode[5][i] = xLd6[i].Average();
                    xLdMode[6][i] = xLd7[i].Average();
                }

                // Calculate means for xLqMode
                for (int i = 0; i < xLq1.Count; ++i)
                {
                    xLqMode[0][i] = xLq1[i].Average();
                    xLqMode[1][i] = xLq2[i].Average();
                    xLqMode[2][i] = xLq3[i].Average();
                    xLqMode[3][i] = xLq4[i].Average();
                    xLqMode[4][i] = xLq5[i].Average();
                    xLqMode[5][i] = xLq6[i].Average();
                    xLqMode[6][i] = xLq6[i].Average();
                }
            }

            else if (type == 2)//xPsimOne[i].Average()
            {
                for (int i = 0; i < xPsimOne.Count; ++i)
                {
                    xPsiMode[0][i] = xPsimOne[i].Max();
                    xPsiMode[1][i] = xPsimTwo[i].Max();
                    xPsiMode[2][i] = xPsimThree[i].Max();
                    xPsiMode[3][i] = xPsimFour[i].Max();
                    xPsiMode[4][i] = xPsimFive[i].Max();
                    xPsiMode[5][i] = xPsimSix[i].Max();
                    xPsiMode[6][i] = xPsimSeven[i].Max();
                }

                // Calculate means for xLdMode
                for (int i = 0; i < xLd1.Count; ++i)
                {
                    xLdMode[0][i] = xLd1[i].Max();
                    xLdMode[1][i] = xLd2[i].Max();
                    xLdMode[2][i] = xLd3[i].Max();
                    xLdMode[3][i] = xLd4[i].Max();
                    xLdMode[4][i] = xLd5[i].Max();
                    xLdMode[5][i] = xLd6[i].Max();
                    xLdMode[6][i] = xLd7[i].Max();
                }

                // Calculate means for xLqMode
                for (int i = 0; i < xLq1.Count; ++i)
                {
                    xLqMode[0][i] = xLq1[i].Max();
                    xLqMode[1][i] = xLq2[i].Max();
                    xLqMode[2][i] = xLq3[i].Max();
                    xLqMode[3][i] = xLq4[i].Max();
                    xLqMode[4][i] = xLq5[i].Max();
                    xLqMode[5][i] = xLq6[i].Max();
                    xLqMode[6][i] = xLq6[i].Max();
                }
            }

            else
            {
                for (int i = 0; i < xPsimOne.Count; ++i)
                {
                    xPsiMode[0][i] = xPsimOne[i][xPsimOne[i].Count - 1];
                    xPsiMode[1][i] = xPsimTwo[i][xPsimTwo[i].Count - 1];
                    xPsiMode[2][i] = xPsimThree[i][xPsimThree[i].Count - 1];
                    xPsiMode[3][i] = xPsimFour[i][xPsimFour[i].Count - 1];
                    xPsiMode[4][i] = xPsimFive[i][xPsimFive[i].Count - 1];
                    xPsiMode[5][i] = xPsimSix[i][xPsimSix[i].Count - 1];
                    xPsiMode[6][i] = xPsimSeven[i][xPsimSeven[i].Count - 1];
                }

                // Fill xLdMode with the last column of respective xLd arrays
                for (int i = 0; i < xLd1.Count; ++i)
                {
                    xLdMode[0][i] = xLd1[i][xLd1[i].Count - 1];
                    xLdMode[1][i] = xLd2[i][xLd2[i].Count - 1];
                    xLdMode[2][i] = xLd3[i][xLd3[i].Count - 1];
                    xLdMode[3][i] = xLd4[i][xLd4[i].Count - 1];
                    xLdMode[4][i] = xLd5[i][xLd5[i].Count - 1];
                    xLdMode[5][i] = xLd6[i][xLd6[i].Count - 1];
                    xLdMode[6][i] = xLd7[i][xLd7[i].Count - 1];
                }



                // Fill xLqMode with the last column of respective xLq arrays
                for (int i = 0; i < xLq1.Count; ++i)
                {
                    xLqMode[0][i] = xLq1[i][xLq1[i].Count - 1];
                    xLqMode[1][i] = xLq2[i][xLq2[i].Count - 1];
                    xLqMode[2][i] = xLq3[i][xLq3[i].Count - 1];
                    xLqMode[3][i] = xLq4[i][xLq4[i].Count - 1];
                    xLqMode[4][i] = xLq5[i][xLq5[i].Count - 1];
                    xLqMode[5][i] = xLq6[i][xLq6[i].Count - 1];
                    xLqMode[6][i] = xLq7[i][xLq7[i].Count - 1];
                }
            }

            List<List<double>> PsiPm = new List<List<double>>(Enumerable.Range(0, 6).Select(x => new List<double>(new double[xPsiMode[0].Count + 1])));

            // Calculate Psi arrays for all 6 rows of the digital twin table
            PsiPm[0] = new List<double> { (Math.Pow(xPsiMode[0][0], 2) / xPsiMode[0][1]) }.Concat(xPsiMode[0]).ToList();
            PsiPm[1] = new List<double> { (Math.Pow(xPsiMode[1][0], 2) / xPsiMode[1][1]) }.Concat(xPsiMode[1]).ToList();
            PsiPm[2] = new List<double> { (Math.Pow(xPsiMode[2][0], 2) / xPsiMode[2][1]) }.Concat(xPsiMode[2]).ToList();
            PsiPm[3] = new List<double> { (Math.Pow(xPsiMode[3][0], 2) / xPsiMode[3][1]) }.Concat(xPsiMode[3]).ToList();
            PsiPm[4] = new List<double> { (Math.Pow(xPsiMode[4][0], 2) / xPsiMode[4][1]) }.Concat(xPsiMode[4]).ToList();
            PsiPm[5] = new List<double> { (Math.Pow(xPsiMode[6][0], 2) / xPsiMode[6][1]) }.Concat(xPsiMode[6]).ToList();

            //std::vector<double> PsiArr(PsiPm[DTnum].size());            //std::vector<double> PsiArr(PsiPm[DTnum].size());

            PLmat["PsiArr"] = new List<List<List<double>>>(Enumerable.Range(0, PsiPm[DTnum].Count).Select(x => new List<List<double>> { new List<double>(new double[1]) }).ToList());
            //PLmat["PsiArr"] = new List<List<double>>(Enumerable.Range(0, PsiPm[DTnum].Count).Select(x => new List<double>(new double[1])));
            //PLmat["PsiArr"] = new OneDList(Enumerable.Range(0, PsiPm[DTnum].Count).Select(x => new List<double>(new double[1])));

            for (int i = 0; i < PsiPm[DTnum].Count; ++i)
            {
                PLmat["PsiArr"][i][0][0] = PsiPm[DTnum][i] * (1 + AdjPM);
            }

            //string testString = PLmat["PsiArr"].Count.ToString() + "x" + PLmat["PsiArr"][0].Count.ToString();
            //_mainWindow.myTextbox.Text = testString;

            List<List<double>> LUTLd = new List<List<double>>(Enumerable.Range(0, 6).Select(x => new List<double>(new double[xLdMode[0].Count + 1])));
            List<List<double>> LUTLq = new List<List<double>>(Enumerable.Range(0, 6).Select(x => new List<double>(new double[xLqMode[0].Count + 1])));


            LUTLd[0] = new List<double> { (Math.Pow(xLdMode[0][0], 2) / xLdMode[0][1]) }.Concat(xLdMode[0]).ToList();
            LUTLd[1] = new List<double> { (Math.Pow(xLdMode[1][0], 2) / xLdMode[1][1]) }.Concat(xLdMode[1]).ToList();
            LUTLd[2] = new List<double> { (Math.Pow(xLdMode[2][0], 2) / xLdMode[2][1]) }.Concat(xLdMode[2]).ToList();
            LUTLd[3] = new List<double> { (Math.Pow(xLdMode[3][0], 2) / xLdMode[3][1]) }.Concat(xLdMode[3]).ToList();
            LUTLd[4] = new List<double> { (Math.Pow(xLdMode[4][0], 2) / xLdMode[4][1]) }.Concat(xLdMode[4]).ToList();
            LUTLd[5] = new List<double> { (Math.Pow(xLdMode[6][0], 2) / xLdMode[6][1]) }.Concat(xLdMode[6]).ToList();

            LUTLq[0] = new List<double> { (Math.Pow(xLqMode[0][0], 2) / xLqMode[0][1]) }.Concat(xLqMode[0]).ToList();
            LUTLq[1] = new List<double> { (Math.Pow(xLqMode[1][0], 2) / xLqMode[1][1]) }.Concat(xLqMode[1]).ToList();
            LUTLq[2] = new List<double> { (Math.Pow(xLqMode[2][0], 2) / xLqMode[2][1]) }.Concat(xLqMode[2]).ToList();
            LUTLq[3] = new List<double> { (Math.Pow(xLqMode[3][0], 2) / xLqMode[3][1]) }.Concat(xLqMode[3]).ToList();
            LUTLq[4] = new List<double> { (Math.Pow(xLqMode[4][0], 2) / xLqMode[4][1]) }.Concat(xLqMode[4]).ToList();
            LUTLq[5] = new List<double> { (Math.Pow(xLqMode[6][0], 2) / xLqMode[6][1]) }.Concat(xLqMode[6]).ToList();


            PLmat["LdArr"] = new List<List<List<double>>>(Enumerable.Range(0, LUTLd[DTnum].Count).Select(x => new List<List<double>> { new List<double>(new double[1]) }).ToList());
            PLmat["LqArr"] = new List<List<List<double>>>(Enumerable.Range(0, LUTLq[DTnum].Count).Select(x => new List<List<double>> { new List<double>(new double[1]) }).ToList());


            //PLmat["LdArr"] = new List<List<double>>(Enumerable.Range(0, LUTLd[DTnum].Count).Select(x => new List<double>(new double[1])));
            //PLmat["LqArr"] = new List<List<double>>(Enumerable.Range(0, LUTLq[DTnum].Count).Select(x => new List<double>(new double[1])));

            for (int i = 0; i < LUTLd[DTnum].Count; ++i)
            {
                PLmat["LdArr"][i][0][0] = LUTLd[DTnum][i] * (1 + AdjLd);
                PLmat["LqArr"][i][0][0] = LUTLq[DTnum][i] * (1 + AdjLq);

            }

            //# Psi_PM sheet
            //    PsiPm = np.zeros((6, xPsiMode.shape[1] + 1))  # Initialize PsiPm array
            //
            //    # Calculate Psi arrays for all 6 rows of the digital twin table
            //    PsiPm[0, :] = np.concatenate([[xPsiMode[0][0] * *2 / xPsiMode[0][1]], xPsiMode[0, :]])  # Col H
            //    PsiPm[1, :] = np.concatenate([[xPsiMode[1][0] * *2 / xPsiMode[1][1]], xPsiMode[1, :]])  # Col I
            //    PsiPm[2, :] = np.concatenate([[xPsiMode[2][0] * *2 / xPsiMode[2][1]], xPsiMode[2, :]])  # Col J
            //    PsiPm[3, :] = np.concatenate([[xPsiMode[3][0] * *2 / xPsiMode[3][1]], xPsiMode[3, :]])  # Col K
            //    PsiPm[4, :] = np.concatenate([[xPsiMode[4][0] * *2 / xPsiMode[4][1]], xPsiMode[4, :]])  # Col L
            //    PsiPm[5, :] = np.concatenate([[xPsiMode[6][0] * *2 / xPsiMode[6][1]], xPsiMode[6, :]])  # Col M
            //
            //    # The selected row is adjusted and returned
            //    # print("PsiPm[DTnum, :]",PsiPm[DTnum, :])
            //    PLmat['PsiArr'] = PsiPm[DTnum, :]* (1 + AdjPM)  # Col C
            //
            //    # LUT(Ld,Lq) sheet
            //    LUTLd = np.zeros((6, xLdMode.shape[1] + 1))  # Initialize LUTLd array
            //    LUTLq = np.zeros((6, xLqMode.shape[1] + 1))  # Initialize LUTLq array
            //
            //    # Calculate Ld arrays for all 6 rows of the digital twin table
            //    LUTLd[0, :] = np.concatenate([[xLdMode[0][0] * *2 / xLdMode[0][1]], xLdMode[0, :]])  # Col D
            //    LUTLd[1, :] = np.concatenate([[xLdMode[1][0] * *2 / xLdMode[1][1]], xLdMode[1, :]])  # Col E
            //    LUTLd[2, :] = np.concatenate([[xLdMode[2][0] * *2 / xLdMode[2][1]], xLdMode[2, :]])  # Col F
            //    LUTLd[3, :] = np.concatenate([[xLdMode[3][0] * *2 / xLdMode[3][1]], xLdMode[3, :]])  # Col G
            //    LUTLd[4, :] = np.concatenate([[xLdMode[4][0] * *2 / xLdMode[4][1]], xLdMode[4, :]])  # Col H
            //    LUTLd[5, :] = np.concatenate([[xLdMode[6][0] * *2 / xLdMode[6][1]], xLdMode[6, :]])  # Col I
            //
            //    # Calculate Lq arrays for all 6 rows of the digital twin table
            //    LUTLq[0, :] = np.concatenate([[xLqMode[0][0] * *2 / xLqMode[0, 1]], xLqMode[0, :]])  # Col O
            //    LUTLq[1, :] = np.concatenate([[xLqMode[1][0] * *2 / xLqMode[1, 1]], xLqMode[1, :]])  # Col P
            //    LUTLq[2, :] = np.concatenate([[xLqMode[2][0] * *2 / xLqMode[2, 1]], xLqMode[2, :]])  # Col Q
            //    LUTLq[3, :] = np.concatenate([[xLqMode[3][0] * *2 / xLqMode[3, 1]], xLqMode[3, :]])  # Col R
            //    LUTLq[4, :] = np.concatenate([[xLqMode[4][0] * *2 / xLqMode[4, 1]], xLqMode[4, :]])  # Col S
            //    LUTLq[5, :] = np.concatenate([[xLqMode[6][0] * *2 / xLqMode[6, 1]], xLqMode[6, :]])  # Col T
            //
            //
            //    # Adjust LdArr and LqArr and return
            //    PLmat['LdArr'] = LUTLd[DTnum, :]* (1 + AdjLd)  # Col J
            //    PLmat['LqArr'] = LUTLq[DTnum, :]* (1 + AdjLq)  # Col U
            //



            //PLmat["DTKk0"] = std::vector<std::vector<std::vector<double>>>(DTKk0.size(), std::vector<std::vector<double>>(DTKk0[0].size(), std::vector<double>(DTKk0[0][0].size(), 0.0)));

            PLmat["DTKk0"] = new List<List<List<double>>>(Enumerable.Range(0, DTKk0.Count).Select(x => new List<List<double>>(Enumerable.Range(0, DTKk0[0].Count).Select(y => new List<double>(new double[DTKk0[0][0].Count])).ToList())).ToList());

            //define a 3D vector obj["DTKk0"] like the 2D vector above
            for (int i = 0; i < DTKk0.Count; i++) //8x27x27
                for (int j = 0; j < DTKk0[0].Count; j++)
                    for (int k = 0; k < DTKk0[0][0].Count; k++)
                        PLmat["DTKk0"][i][j][k] = DTKk0[i][j][k];

            //PLmat["DTKk0"] = DTKk0;
            PLmat["Hdmat"] = new List<List<List<double>>>(Enumerable.Range(0, Hdmat.Count).Select(x => new List<List<double>>(Enumerable.Range(0, Hdmat[0].Count).Select(y => new List<double>(new double[1])).ToList())).ToList());
            PLmat["Hqmat"] = new List<List<List<double>>>(Enumerable.Range(0, Hqmat.Count).Select(x => new List<List<double>>(Enumerable.Range(0, Hqmat[0].Count).Select(y => new List<double>(new double[1])).ToList())).ToList());
            PLmat["Ld"] = new List<List<List<double>>>(Enumerable.Range(0, LUTLd.Count).Select(x => new List<List<double>>(Enumerable.Range(0, LUTLd[0].Count).Select(y => new List<double>(new double[1])).ToList())).ToList());
            PLmat["Lq"] = new List<List<List<double>>>(Enumerable.Range(0, LUTLq.Count).Select(x => new List<List<double>>(Enumerable.Range(0, LUTLq[0].Count).Select(y => new List<double>(new double[1])).ToList())).ToList());
            PLmat["Psi"] = new List<List<List<double>>>(Enumerable.Range(0, PsiPm.Count).Select(x => new List<List<double>>(Enumerable.Range(0, PsiPm[0].Count).Select(y => new List<double>(new double[1])).ToList())).ToList());


            for (int i = 0; i < Hdmat.Count; i++) //6x27
                for (int j = 0; j < Hdmat[0].Count; j++)
                {
                    //cout << "ixj" << i << j;
                    PLmat["Hdmat"][i][j][0] = Hdmat[i][j];
                    PLmat["Hqmat"][i][j][0] = Hqmat[i][j];

                }

            for (int i = 0; i < LUTLd.Count; i++) //6x28
                for (int j = 0; j < LUTLd[0].Count; j++)
                {
                    PLmat["Ld"][i][j][0] = LUTLd[i][j];
                    PLmat["Lq"][i][j][0] = LUTLq[i][j];
                    PLmat["Psi"][i][j][0] = PsiPm[i][j];
                }
            
            PLmat["k0mat"] = new List<List<List<double>>>(Enumerable.Range(0, 6).Select(x => new List<List<double>>(Enumerable.Range(0, oneK0.Count).Select(y => new List<double>(new double[oneK0[0].Count])).ToList())).ToList());
            for (int i = 0; i < oneK0.Count; i++)
                for (int j = 0; j < oneK0[0].Count; j++)
                {
                PLmat["k0mat"][0][i][j] = oneK0[i][j];
                PLmat["k0mat"][1][i][j] = twoK0[i][j];
                PLmat["k0mat"][2][i][j] = threeK0[i][j];
                PLmat["k0mat"][3][i][j] = fourK0[i][j];
                PLmat["k0mat"][4][i][j] = fiveK0[i][j];
                PLmat["k0mat"][5][i][j] = sevenK0[i][j];
            }


            return PLmat;
        }



        public void startupFcn(double pole, double flux, double Nominal_d_axis_Inductance, double Nominal_q_axis_Inductance, double Resistance, double Voltage, double cur)
        {

            // Initialize motor parameters
            mot = new Dictionary<string, double>
        {
            { "pole", pole },
            { "flux", flux },
            { "res", Resistance },
            { "cfe", 0.24 },
            { "cstr", 0.6 },
            { "Cstrunit", 0.000000001 },
            { "cPfric", 3 },
            { "cPwind", 2 },
            { "freq", 6500 },
            { "vol", Voltage },
            { "cur", cur },
            { "DCV", 350 },
            { "DCC", 480 },
            { "unit", 0.000001 },
            { "ld", Nominal_d_axis_Inductance },
            { "lq", Nominal_q_axis_Inductance }
        };

            EffMaxSpeed = 10000;
            mot["d"] = (mot["lq"] - mot["ld"]) * mot["unit"];
            mot["e"] = mot["lq"] / mot["ld"];
            mot["If"] = mot["flux"] / (mot["ld"] * mot["unit"]);

            double vm = 350;
            //double[,] tabValues = new double[,]
            //{
            //{ 0.1370, 0.2110, 222, 323, 0.085, 0, 1.2800 },
            //{ 0.1170, 0.1820, 230, 333, 0.092, 0, 1.18 },
            //{ 0.1614, 0.184, 289.3, 382.8, 0.1262, 0, 1 },
            //{ 0.106, 0.1610, 276, 403, 0.1034, 0, 1.0600 },
            //{ 0.1510, 0.2250, 250, 366, 0.0907, 0, 1.3900 },
            //{ 0.141, 0.229, 212.7, 311.2, 0.07296, 0, 1.4100 }
            //};

            List<List<double>> tabValues = new List<List<double>>()
            {
                new List<double> { 0.1370, 0.2110, 222, 323, 0.085, 0, 1.2800 },
                new List<double> { 0.1170, 0.1820, 230, 333, 0.092, 0, 1.18 },
                new List<double> { 0.1614, 0.184, 289.3, 382.8, 0.1262, 0, 1 },
                new List<double> { 0.106, 0.1610, 276, 403, 0.1034, 0, 1.0600 },
                new List<double> { 0.1510, 0.2250, 250, 366, 0.0907, 0, 1.3900 },
                new List<double> { 0.141, 0.229, 212.7, 311.2, 0.07296, 0, 1.4100 }
            };


            double[] ldopt = new double[6];
            double[] lqopt = { 87, 88, 120, 100, 90, 76 };
            double[] pmdopt = { 18, 13, 0, 20, 11, 0 };

            //object[] lqoptObjectArray = lqopt.Cast<object>().ToArray();
            //string arrayString = Convert1DArrayToString(lqoptObjectArray);
            //_mainWindow.myTextbox.Text = arrayString;

            List<List<double>>DTtab =new List<List<double>>(); 
            DTtab = DigiTwin(mot["pole"], vm, tabValues);

            //var ploeString= mot["pole"].ToString();
            //_mainWindow.myTextbox.Text = ploeString;

            //var DTtabstring = Convert2DListToString(DTtab);
            //_mainWindow.myTextbox.Text = DTtabstring;

            // Read the specified sheet
            DataTable a12efftab = ReadExcelSheet("InputTableFile.xlsx", "A12 Efficiency");

            // Fill missing (NaN) values with 0, similar to fillna(0) in pandas
            foreach (DataRow row in a12efftab.Rows)
            {
                for (int i = 0; i < a12efftab.Columns.Count; i++)
                {
                    if (row[i] == DBNull.Value)
                    {
                        row[i] = 0;
                    }
                }
            }

            //int rows = a12efftab.Rows.Count;
            //int columns = a12efftab.Columns.Count;

            // Create a 2D array with the same dimensions as the DataTable
            //object[,] array2D = new object[rows, columns];

            object[] prpm = new object[a12efftab.Rows.Count];
            object[] posId = new object[a12efftab.Rows.Count];
            object[] posIq = new object[a12efftab.Rows.Count];
            object[] nrpm = new object[a12efftab.Rows.Count];
            object[] negId = new object[a12efftab.Rows.Count];
            object[] negIq = new object[a12efftab.Rows.Count];
            object[] Te = new object[a12efftab.Rows.Count];
            object[] I = new object[a12efftab.Rows.Count];

            // Copy the DataTable data into the array
            for (int i = 1; i < a12efftab.Rows.Count; i++) 
            {
                prpm[i-1] = a12efftab.Rows[i][0];
                posId[i-1] = a12efftab.Rows[i][1];
                posIq[i - 1] = a12efftab.Rows[i][2];
                nrpm[i - 1] = a12efftab.Rows[i][10];
                negId[i - 1] = a12efftab.Rows[i][11];
                negIq[i - 1] = a12efftab.Rows[i][12]; 
                Te[i - 1] = a12efftab.Rows[i][13];;
                I[i - 1] = a12efftab.Rows[i][15];
            }

            var a12 = new Dictionary<string, object>
            {
                { "prpm", prpm },
                { "posId", posId },
                { "posIq", posIq },
                { "nrpm", nrpm },
                { "negId", negId },
                { "negIq", negIq },
                { "Te", Te },
                { "I", I }
            };



            //string dictionaryString = ConvertDictionaryKeyToString(a12, "posId");
            //_mainWindow.myTextbox.Text = dictionaryString;

            int type = 2;
            int hdcalf = 1;
            int hqcalf = 1;
            int _const = 0;
            int tnflag = 1;
            int peflag = 1;
            int mode = 2;
            int DTnum = 2;

            DataTable _Psid = ReadExcelSheet("InputTableFile.xlsx", "Psi_d");
            DataTable _Psiq = ReadExcelSheet("InputTableFile.xlsx", "Psi_q");

            //var Psidstring = DataTableToString(_Psid);
            //_mainWindow.myTextbox.Text = Psidstring;

            List<List<double>> Psid = new List<List<double>>(Enumerable.Range(0, _Psid.Rows.Count - 1).Select(x => new List<double>(new double[_Psid.Columns.Count - 2])));

            //double[,] Psid = new double[_Psid.Rows.Count-1, _Psid.Columns.Count-2];


            for (int i = 1; i < _Psid.Rows.Count; i++)  //Skip one row
                for (int j = 2; j < _Psid.Columns.Count; j++)   //Skip one column
                {
                    Psid[i - 1][j - 2] = (double)_Psid.Rows[i][j];
                }

            //double[,] Psiq = new double[_Psiq.Rows.Count - 1, _Psiq.Columns.Count - 2];

            List<List<double>> Psiq = new List<List<double>>(Enumerable.Range(0, _Psiq.Rows.Count - 1).Select(x => new List<double>(new double[_Psiq.Columns.Count - 2])));

            for (int i = 1; i < _Psiq.Rows.Count; i++)  //Skip one row
                for (int j = 2; j < _Psiq.Columns.Count; j++)   //Skip one column
                {
                    Psiq[i - 1][j - 2] = (double)_Psiq.Rows[i][j];
                }

            //var Psiqstring = Convert2DArrayToString(Psiq);
            //_mainWindow.myTextbox.Text = Psiqstring;

            //double[,] Omega = new double[Psiq.GetLength(0), Psiq.GetLength(1)];
            List<List<double>> Omega = new List<List<double>>(Enumerable.Range(0, Psiq.Count).Select(x => new List<double>(new double[Psiq[0].Count])));

            //_mainWindow.myTextbox.Text = Psiq.Count.ToString();
            //_mainWindow.myTextbox.Text += "," + Psiq[0].Count.ToString();

            for (int i = 0; i < Omega.Count; i++)  //Skip one row
                for (int j = 0; j < Omega[0].Count; j++)   //Skip one column
                {
                    Omega[i][j] = Math.Pow(Psid[i][j], 2) + Math.Pow(Psiq[i][j], 2);
                }




            Dictionary<string, List<List<List<double>>>> PLmat = PsiLdLq(mot["pole"], DTtab, hdcalf, hqcalf, type, DTnum, pmdopt[DTnum], ldopt[DTnum], lqopt[DTnum], Psid, Psiq, Omega);

            DataTable tunTable = ReadExcelSheet("InputTableFile.xlsx", "tuning");
             
            PLmat["copy"] = new List<List<List<double>>>(Enumerable.Range(0, tunTable.Rows.Count).Select(x => new List<List<double>> { new List<double>(new double[1]) }).ToList());  //1D
            PLmat["tuning"] = new List<List<List<double>>>(Enumerable.Range(0, tunTable.Rows.Count).Select(x => new List<List<double>> { new List<double>(new double[1]) }).ToList());  //1D

            for (int i = 0; i < PLmat["copy"].Count; i++)  //Skip one row
                PLmat["copy"][i][0][0] = (double)tunTable.Rows[i][1];

            //var testString = Convert3DListToString(PLmat["copy"]);
            //_mainWindow.myTextbox.Text = testString;

            //_mainWindow.myTextbox.Text = PLmat["copy"].Count.ToString();


            int tun =1;
            if (tun==1)
            {
                for (int i = 0; i < PLmat["tuning"].Count; i++)  //Skip one row
                    PLmat["tuning"][i][0][0] = (double)tunTable.Rows[i][0];
            }
            else
            {
                for (int i = 0; i < PLmat["tuning"].Count; i++)  //Skip one row
                    PLmat["tuning"][i][0][0] = 1;
            }

            double S1K = 1000;
            double S5K = 5000;
            double S10K = 10000;
            double S15K = 15000;

            ipmclass ipm = new ipmclass();

            ipm.initial = 0.0001;
            ipm.increment = 200;

            //measuring the execution time of ipm init function
            Stopwatch stopwatch = new Stopwatch(); // Create a new Stopwatch instance
            stopwatch.Start(); // Start measuring time
            
            ipm.init(mot, DTtab, DTnum - 1, PLmat, tnflag, peflag, _const);
            //foreach (var value in ipm.Tn)
            //{
            //    Console.WriteLine(value);
            //}

            stopwatch.Stop(); // Stop measuring time
            Console.WriteLine($"Execution Time: {stopwatch.ElapsedMilliseconds} ms");

            /////////////////////////////wltc///////////

            double lowf = 1.2;
            double midf = 1.0;
            double highf = 1.0;
            double exhif = 1.0;


            DataTable _tbl = ReadExcelSheet("InputTableFile.xlsx", "Cruising Pattern");
            var ip = new Dictionary<string, dynamic>
            {
                { "gravity", 9.81 },
                { "fr", 0.008 },
                { "Vwind", 0.0 },
                { "Cd", 0.25 },
                { "Af", 1.5 },
                { "p", 1.225 },
                { "rdgrd", 0.0 },
                { "gdr", 10.0 },
                { "regen_ratio", 0.8 },
                { "regen_limit", 16.0 },
                {"GrossWt",2000.0},
                {"J",272.3586182994925},
                {"zone",2 },
               // { "VDswitch",1}
            };
            string columnName = "";
            switch ((int)ip["zone"])
            {
                case 1:
                    columnName = "NEDC";
                    break;
                case 2:
                    columnName = "WLTC";
                    break;
                case 3:
                    columnName = "CLTC";
                    break;
                case 4:
                    columnName = "FTP75";
                    break;
                default:
                    throw new Exception("Invalid zone value.");
            }
            List<double> crsPtn = new List<double>();
            foreach (DataRow row in _tbl.Rows)
            {
                if (row[columnName] != DBNull.Value)
                {
                    crsPtn.Add(Convert.ToDouble(row[columnName]));
                }
            }
            ip["crsPtn"] = crsPtn;


            wltcclass dp = new wltcclass();
            dp.sec = new List<double>();//1800
            var kmh= new List<double>(); ;//1800
            DataTable _wltctable = ReadExcelSheet("InputTableFile.xlsx", "WLTC-Class3");
            
            string wltcColName = "S";
            foreach (DataRow row in _wltctable.Rows)
            {
                if (row[wltcColName] != DBNull.Value)
                {
                    dp.sec.Add(Convert.ToDouble(row[wltcColName]));
                }
            }
            string wltcColName1 = "km_h";
            //dp.num = len(kmh)
            foreach (DataRow row in _wltctable.Rows)
            {
                if (row[wltcColName1] != DBNull.Value)
                {
                    kmh.Add(Convert.ToDouble(row[wltcColName1]));
                }
            }
            dp.num = kmh.Count;
            double wlow = 590.0;
            double wmiddle = 1023.0;
            double whigh = 1478.0;
            double sx = 0.02;
            bool Vd = true;
            List<double> Vx = new List<double>((int)dp.num);
            Console.WriteLine("dp.sec length: " +Vx.Count);
            for (int i = 0; i < dp.num; i++)
            {
                Vx.Add(0.0); // Initialize with default value (0.0)
            }

            for (int idx = 0; idx < dp.num; idx++)
            {
                if (Vd)
                {
                    ip["VDswitch"] = 1;
                    if (idx < wlow)
                    {
                       // Console.WriteLine(crsPtn[idx]);
                        Vx[idx] = crsPtn[idx] * lowf / (1 - sx);
                    }
                    else if (idx < wmiddle)
                    {
                        Vx[idx] = crsPtn[idx] * midf / (1 - sx);
                    }
                    else if (idx < whigh)
                    {
                        Vx[idx] = crsPtn[idx] * highf / (1 - sx);
                    }
                    else
                    {
                        Vx[idx] = crsPtn[idx] * exhif / (1 - sx);
                    }
                }
                else
                {
                    ip["VDswitch"] = 0;
                    if (idx < wlow)
                    {
                        // Console.WriteLine(crsPtn[idx]);
                        Vx[idx] = crsPtn[idx] * lowf  ;
                    }
                    else if (idx < wmiddle)
                    {
                        Vx[idx] = crsPtn[idx] * midf  ;
                    }
                    else if (idx < whigh)
                    {
                        Vx[idx] = crsPtn[idx] * highf ;
                    }
                    else
                    {
                        Vx[idx] = crsPtn[idx] * exhif  ;
                    }

                }

            }
            ip["Vx"] = Vx;
            var btt = new Dictionary<string, double>
            {
                { "charge",80 },
                { "regen_ratio", 0.750011315986961 },
                { "regen_limit",16.825305674317924},
            };
            var gear = new Dictionary<string, double>
            {
                { "gdr",10 },
                { "shaftD", 0.07 },

            };
            var tire = new Dictionary<string, double>
            {
                { "Ma",10 },
                { "TireOutD", 29.680044047914375 },
                { "Tirerw", 0.37693655940851256},
            };

            dp.init(mot["pole"], btt, gear, tire, ip, ipm);


            //  var crsPtnList = (List<double>)ip["crsPtn"]; // Cast the object to List<double>
            //Console.WriteLine("dp.sec length: " + dp.sec.Count);
            //foreach (var value in ipm.Tn)
            //{
            //    Console.WriteLine(value);
            //}

            // List<double> crsPtn = GetCruisingPattern(tbl, (int)ip["zone"]);

            //{
            //    Initial = 0.0001,
            //    Increment = 200
            //};
            //    //    ipm.Init(mot, DTtab, DTnum - 1, PLmat, tnflag, peflag, 0);

            //    //    double[] speeds = { 806.201, 804.691, 796.900, 816.030, 836.279, 807.259 };

            //    //    inv = new Dictionary<string, double>
            //    //{
            //    //    { "tr", 0.09 },
            //    //    { "tf", 0.16 },
            //    //    { "ton", 124.8 },
            //    //    { "von", 0.7 },
            //    //    { "trr", 0.06 }
            //    //};

            //    //    igbt = new Dictionary<string, double>
            //    //{
            //    //    { "modulation", 0.8 },
            //    //    { "A1", 0.0012 },
            //    //    { "A0", 0.8965 },
            //    //    { "C3", 3e-8 },
            //    //    { "C2", -0.000028 },
            //    //    { "C1", 0.058 },
            //    //    { "C0", 1.365 },
            //    //    { "D3", 7.4e-8 },
            //    //    { "D2", -0.000084 },
            //    //    { "D1", 0.09 },
            //    //    { "D0", 4.021 },
            //    //    { "B1", 0.0015 },
            //    //    { "B0", 0.84 },
            //    //    { "E3", -2.9e-8 },
            //    //    { "E2", 0.000025 },
            //    //    { "E1", 0.034 },
            //    //    { "E0", 7.016 }
            //    //};

            //    //    Temp = new Dictionary<string, double>
            //    //{
            //    //    { "A", 0.004 },
            //    //    { "B", 0.0005 },
            //    //    { "C", 0.0566 },
            //    //    { "iniTemp", 65 },
            //    //    { "RoomT", 23 }
            //    //};

            //    //    Flag = new Dictionary<string, int>
            //    //{
            //    //    { "Pcu", 1 },
            //    //    { "Pfe", 2 },
            //    //    { "Pstr", 1 },
            //    //    { "Pf", 1 },
            //    //    { "Pw", 1 },
            //    //    { "Pinv", 2 },
            //    //    { "temp", 1 }
            //    //};

            //    //    tire = new Dictionary<string, double>
            //    //{
            //    //    { "Ma", 10 },
            //    //    { "TireOutD", 27 }
            //    //};

            //    //    vhcl = new Dictionary<string, double>
            //    //{
            //    //    { "crb", 1500 },
            //    //    { "psgw", 80 },
            //    //    { "psgnum", 4 },
            //    //    { "lug", 140 },
            //    //    { "tireno", 4 }
            //    //};

            //    //    gear = new Dictionary<string, double>
            //    //{
            //    //    { "gdr", 10 },
            //    //    { "shaftD", 0.07 }
            //    //};

            //    //    ip = new Dictionary<string, double>
            //    //{
            //    //    { "gravity", 9.81 },
            //    //    { "fr", 0.008 },
            //    //    { "Vwind", 0 },
            //    //    { "Cd", 0.25 },
            //    //    { "Af", 1.5 },
            //    //    { "p", 1.225 },
            //    //    { "rdgrd", 0 },
            //    //    { "gdr", 10 },
            //    //    { "regen_ratio", 0.8 },
            //    //    { "regen_limit", 16 }
            //    //};

            //    //    double sx = 0.02;

            //    //    double lowf = 1.2;
            //    //    double midf = 1.0;
            //    //    double highf = 1.0;
            //    //    double exhif = 1.0;

            //    //    double Tn50 = 50;
            //    //    double Tn100 = 100;
            //    //    double Tn200 = 200;
            //    //    double Tn300 = 300;
            //    //}

            //    //private static DataTable ReadExcelSheet(string fileName, string sheetName, bool withHeader = false)
            //    //{
            //    //    // Implement reading from Excel file here (e.g., using ExcelDataReader or another library).
            //    //    // Returning a DataTable as an example.
            //    //    return new DataTable();
            //    //}
        }
    }
    }
