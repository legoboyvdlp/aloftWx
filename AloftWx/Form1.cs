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

// This is the code for your desktop app.
// Press Ctrl+F5 (or go to Debug > Start Without Debugging) to run your app.

namespace AloftWx
{
    public partial class AloftWxForm : Form
    {
        public AloftWxForm()
        {
            InitializeComponent();
        }

        private void linkLabel1_LinkClicked(object sender, LinkLabelLinkClickedEventArgs e)
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

            if (int.Parse(forecast) < 10)
            {
                forecast = "0" + forecast;
            }

            // string NOAAurl = "http://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/gfs." + cycleStr + "/" + filename;
            string NOAAurl = "https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p50.pl?file=gfs.t" + curCycleStr + "z.pgrb2full.0p50.f" + forecast + "&lev_150_mb=on&lev_175_mb=on&lev_200_mb=on&lev_225_mb=on&lev_250_mb=on&lev_275_mb=on&lev_300_mb=on&lev_325_mb=on&lev_350_mb=on&lev_400_mb=on&lev_450_mb=on&lev_500_mb=on&lev_550_mb=on&lev_600_mb=on&lev_650_mb=on&lev_700_mb=on&lev_750_mb=on&lev_850_mb=on&lev_900_mb=on&var_UGRD=on&var_VGRD=on&leftlon=0&rightlon=360&toplat=90&bottomlat=-90&dir=%2Fgfs." + cycleStr;
            
            urlBox.Text = NOAAurl;

            // call file dialog
            SaveFileDialog SaveFileDialog1 = new SaveFileDialog();
            SaveFileDialog1.Filter = "GRIB2|*.grib2";
            SaveFileDialog1.Title = "Save GRIB2 file";
            SaveFileDialog1.ShowDialog();

            if (SaveFileDialog1.FileName != "")
            {
                DownloadCurrent(NOAAurl, SaveFileDialog1.FileName);
            }

            //Program.LaunchWGrib2();
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

        private void Form1_Load(object sender, EventArgs e)
        {
           
        }
    }
}
