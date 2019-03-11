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

        private void Button1_Click(object sender, EventArgs e)
        {
            string cycleStr = Program.GetDateCycle();
            int curCycle = int.Parse(cycleStr.Substring(Math.Max(0, cycleStr.Length - 2)));
            string forecast = Program.GetForecast(curCycle);
            string filename = Program.GetFilename(cycleStr, forecast);
            
            // HTTP request

            Stream myStream;
            SaveFileDialog saveFileDialog1 = new SaveFileDialog();

            saveFileDialog1.Filter = "All files (*.*)|*.*";
            saveFileDialog1.FilterIndex = 2;
            saveFileDialog1.RestoreDirectory = true;

            if (saveFileDialog1.ShowDialog() == DialogResult.OK)
            {
                if ((myStream = saveFileDialog1.OpenFile()) != null)
                {
                    string NOAAurl = "http://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/gfs." + cycleStr + "/" + filename;
                    Console.WriteLine(NOAAurl);
                    HttpWebRequest WebReq = (HttpWebRequest)WebRequest.Create
                        (string.Format(NOAAurl));

                    WebReq.Method = "GET";

                    HttpWebResponse WebResp = (HttpWebResponse)WebReq.GetResponse();
                    if (WebResp.StatusCode != HttpStatusCode.OK)
                    {
                        Console.WriteLine("Server status: " + WebResp.StatusCode);
                    }
                    Stream Answer = WebResp.GetResponseStream();
                    Answer.CopyTo(myStream);
                    
                    myStream.Close();
                }
            }
        }

        private void Form1_Load(object sender, EventArgs e)
        {
           
        }
    }
}
