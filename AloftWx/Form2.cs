using AloftWx.Properties;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace AloftWx
{
    public partial class Form2 : Form
    {
        public Form2()
        {
            InitializeComponent();
        }

        private void Button1_Click(object sender, EventArgs e)
        {
            FolderBrowserDialog OpenFolderDialog2 = new FolderBrowserDialog();
            DialogResult result = OpenFolderDialog2.ShowDialog();
            if (result == DialogResult.OK && !string.IsNullOrWhiteSpace(OpenFolderDialog2.SelectedPath))
            {
                wgrib2pathbox.Text = OpenFolderDialog2.SelectedPath;
            }
        }

        private void button2_Click(object sender, EventArgs e)
        {
            { Settings.Default.PATHwgrib2 = wgrib2pathbox.Text; }
            this.Close();
        }

        private void udpButton_Click(object sender, EventArgs e)
        {
            if (udpBox.Text != "")
            {
                { Settings.Default.UDPportOut = int.Parse(udpBox.Text); }
            }
        }

        private void udpButtonIn_Click(object sender, EventArgs e)
        {
            if (udpBoxIn.Text != "")
            {
                { Settings.Default.UDPportIn = int.Parse(udpBoxIn.Text); }
            }
        }

        private void Form2_Closed(object sender, FormClosedEventArgs e)
        {
            if (wgrib2pathbox.Text != "")
            {
                { Settings.Default.PATHwgrib2 = wgrib2pathbox.Text; }
            }

            if (udpBox.Text != "")
            {
                { Settings.Default.UDPportOut = int.Parse(udpBox.Text); }
            }

            if (udpBoxIn.Text != "")
            {
                { Settings.Default.UDPportIn = int.Parse(udpBoxIn.Text); }
            }

            Settings.Default.Save();
        }

        private void Form2_Load(object sender, EventArgs e)
        {
            if (wgrib2pathbox.Text != Settings.Default.PATHwgrib2)
            {
                { wgrib2pathbox.Text = Settings.Default.PATHwgrib2; }
            }

            if (udpBox.Text != Settings.Default.UDPportOut.ToString())
            {
                {udpBox.Text = Settings.Default.UDPportOut.ToString(); }
            }

            if (udpBoxIn.Text != Settings.Default.UDPportIn.ToString())
            {
                { udpBoxIn.Text = Settings.Default.UDPportIn.ToString(); }
            }

            this.FormClosed += new FormClosedEventHandler(Form2_Closed);
        }
    }
}
