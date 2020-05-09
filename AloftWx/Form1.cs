using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;
using System.Windows.Forms;
using System.IO;
using System.Net;
using AloftWx.Properties;
using System.Net.Sockets;

namespace AloftWx
{
    public partial class AloftWxForm : Form
    {
        private string savedFileName = "";
        private bool _connectFlag = false;
        private bool _parseFlag = false;

        public AloftWxForm()
        {
            InitializeComponent();
        }

        private void LinkLabel1_LinkClicked(object sender, LinkLabelLinkClickedEventArgs e)
        {
            // Click on the link below to continue learning how to build a desktop app using WinForms!
            System.Diagnostics.Process.Start("http://aka.ms/dotnet-get-started-desktop");

        }

        private void Button1_ClickAsync(object sender, EventArgs e)
        {
            string cycleStr = Program.GetDateCycle();
            int curCycle = int.Parse(cycleStr.Substring(Math.Max(0, cycleStr.Length - 2)));
            string curCycleStr = curCycle.ToString();

            if (curCycle < 10)
            {
                curCycleStr = "0" + curCycleStr;
            }

            string forecast = Program.GetForecast(curCycle);
            string filename = Program.GetFilename(cycleStr, forecast);

            string NOAAurl = "https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p50.pl?file=gfs.t" + curCycleStr + "z.pgrb2full.0p50.f" + forecast + "&lev_100_mb=on&lev_125_mb=on&lev_150_mb=on&lev_175_mb=on&lev_200_mb=on&lev_225_mb=on&lev_250_mb=on&lev_275_mb=on&lev_300_mb=on&lev_325_mb=on&lev_350_mb=on&lev_375_mb=on&lev_400_mb=on&lev_425_mb=on&lev_450_mb=on&lev_475_mb=on&lev_500_mb=on&lev_525_mb=on&lev_550_mb=on&lev_575_mb=on&lev_600_mb=on&lev_625_mb=on&lev_650_mb=on&lev_675_mb=on&lev_700_mb=on&lev_725_mb=on&lev_750_mb=on&lev_775_mb=on&lev_800_mb=on&lev_825_mb=on&lev_850_mb=on&lev_875_mb=on&lev_900_mb=on&lev_925_mb=on&lev_950_mb=on&lev_975_mb=on&lev_1000_mb=on&var_UGRD=on&var_VGRD=on&leftlon=0&rightlon=360&toplat=90&bottomlat=-90&dir=%2Fgfs." + cycleStr;
            
            urlBox.Text = NOAAurl;

            // call file dialog
            SaveFileDialog SaveFileDialog1 = new SaveFileDialog
            {
                Filter = "GRIB2|*.grib2",
                Title = "Save GRIB2 file",
                //SaveFileDialog1.InitialDirectory = Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);
                InitialDirectory = Application.StartupPath,
                RestoreDirectory = true,
                FileName = "gfs.t" + curCycleStr + "z0p50f" + forecast + ".grib2"
            };
            SaveFileDialog1.ShowDialog();

            if (SaveFileDialog1.FileName != "")
            {
                DownloadCurrent(NOAAurl, SaveFileDialog1.FileName);
                savedFileName = SaveFileDialog1.FileName;
            }
        }
        
        private void DownloadCurrent(string url,string filename)
        {
            var wc = new WebClient();
            wc.DownloadProgressChanged += (s, e) =>
            {
                 UpdateBar(e.ProgressPercentage);
            };
            wc.DownloadDataCompleted += (sender, e) => {
                if (e.Error != null)
                {
                    Console.WriteLine("Error when downloading file: ", e.Error);
                }
                File.WriteAllBytes(filename, e.Result);
            };
            wc.DownloadDataAsync(new Uri(url));
        }

        public void UpdateBar(int value)
        {
            progressBar1.Value = value;
        }

        public void UpdateBar2(int value)
        {
            progressBar2.Value = value;
        }

        private void Opengrib2FileToolStripMenuItem_Click(object sender, EventArgs e)
        {
            OpenFileDialog OpenFileDialog = new OpenFileDialog
            {
                Title = "Choose grib2 file"
            };
            OpenFileDialog.ShowDialog();
            
            loadedfileBox.Text = OpenFileDialog.FileName;
        }

        private void Button2_Click(object sender, EventArgs e)
        {
            if (savedFileName == "")
            {
                MessageBox.Show("No file downloaded", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                return;
            }
            loadedfileBox.Text = savedFileName;
        }
        
        private void QuitToolStripMenuItem_Click(object sender, EventArgs e)
        {
            this.Close();
        }

        private void Button1_Click(object sender, EventArgs e)
        {
            ParseGribBtn();
        }

        private bool ParseGribBtn()
        {
            if (_parseFlag == false)
            {
                int progress1 = Program.LaunchWGrib2(loadedfileBox.Text);
                UpdateBar2(progress1);
                _parseFlag = true;
            }
            return _parseFlag;
        }

        private void Button3_Click(object sender, EventArgs e)
        {
            OpenFileDialog OpenFileDialog = new OpenFileDialog
            {
                Title = "Choose grib2 file"
            };
            OpenFileDialog.ShowDialog();

            loadedfileBox.Text = OpenFileDialog.FileName;
        }

        private void OptionsToolStripMenuItem_Click_1(object sender, EventArgs e)
        {
            Form2 f2 = new Form2();
            f2.ShowDialog(); // Shows Form2
        }
        
        private void ConnectToFGToolStripMenuItem_Click(object sender, EventArgs e)
        {
            StartOutConnection();
        }

        private bool StartOutConnection()
        {
            if (_connectFlag == true) {
                return _connectFlag = true;
            }
            statusLabel.Text = "Connection Open";
            Program.ReceiveData(Settings.Default.UDPportOut);
            return _connectFlag = true;
        }
        
    }
}
