using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.ConstrainedExecution;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;

namespace Genetic_WPF
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {

        public MainWindow()
        {

            // Register the encoding provider
            System.Text.Encoding.RegisterProvider(System.Text.CodePagesEncodingProvider.Instance);

            GeneticAlgorithm.mot["hello"] = 1.3276768;
            InitializeComponent();
            //myTextbox.Text = GeneticAlgorithm.mot["hello"].ToString("F2");

            //(double pole, double flux, double Nominal_d_axis_Inductance, double Nominal_q_axis_Inductance, double Resistance, double Voltage, double cur)

            GeneticAlgorithm app = new GeneticAlgorithm(this);

            //def startupFcn(pole, flux, Nominal_d_axis_Inductance, Nominal_q_axis_Inductance, Resistance, Voltage, cur):
            //app.startupFcn(6, 0.1262, 289.3, 382.8, 0.0125, 350, 680.9480736639688);
            app.startupFcn(6, 0.1262, 289.3, 382.8, 0.0125, 358.8309291392894, 680.9480736639688);

            //myTextbox.Text="a12 prpm value";
        }

        // Property to access myTextbox
        public TextBox MyTextbox
        {
            get { return myTextbox; }
        }
    }
}
