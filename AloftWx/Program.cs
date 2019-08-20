using AloftWx.Properties;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Net;
using System.Net.Sockets;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace AloftWx
{
    static class Program
    {
        private static string PATHwgrib2;

        public static double alt = -9999;
        public static double lon = -9999;
        public static double lat = -9999;
        public static double[] windSpd = new double[9];
        public static double[] windDir = new double[9];

        private static bool _willParseflag = false;
        private static double latRounded;
        private static double lonRounded;

        [STAThread]
        static void Main()
        {
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            Application.Run(new AloftWxForm());
        }

        public static string GetDateCycle()
        {
            // declare variables
            long currentYear = long.Parse(DateTime.UtcNow.Year.ToString());
            long currentMonth = long.Parse(DateTime.UtcNow.Month.ToString());
            long currentDay = long.Parse(DateTime.UtcNow.Day.ToString());
            long currentHour = long.Parse(DateTime.UtcNow.Hour.ToString());
            long curCycle = 0;

            // declare arrays
            int[] cycles = new int[] { 0, 6, 12, 18 };

            // adjust hour by removing 5
            // a cycle is published between 4 and 5 hours after the cycle time. 
            // cycle gfs.2018121806 was published at 10:54 UTC (4 hours 54 min)
            currentHour = currentHour - 5;

            // if the adjusted cycle is less than 0, add 24 hours and reduce day by 1
            if (currentHour < 0)
            {
                currentHour += 24;
                currentDay -= 1;
            }

            // get the correct cycle
            for (var i = 3; i >= 0; i--)
            {
                if (currentHour >= cycles[i])
                {
                    curCycle = cycles[i];
                    break;
                }
            }

            // if the day is less than 10 append an extra 0
            string curDayStr = currentDay.ToString();
            if (currentDay < 10)
            {
                curDayStr = 0 + curDayStr;
            }

            // if the cycle is 0, add an extra 0
            string curCycleStr = curCycle.ToString();
            if (curCycle < 10)
            {
                curCycleStr = 0 + curCycleStr;
            }

            // if the month is less than 10, add an extra 0
            string currentMonthStr = currentMonth.ToString();
            if (currentMonth < 10) { currentMonthStr = 0 + currentMonth.ToString(); }
            return currentYear.ToString() + currentMonthStr + curDayStr + curCycleStr;
        }

        public static string GetForecast(int cycle)
        {
            // declare variable
            long currentHour = long.Parse(DateTime.UtcNow.Hour.ToString());
            long currentDay = long.Parse(DateTime.UtcNow.Day.ToString());
            int forecast = 0;

            // declare array
            int[] forecasts = new int[] { 6, 9, 12, 15, 18, 21, 24 };

            // adjust hour by removing 5
            currentHour = currentHour - 5;

            // if the adjusted hour is less than 0, add 24 hours and reduce day by 1
            if (currentHour < 0)
            {
                currentHour += 24;
                currentDay -= 1;
            }

            // get the correct forecast
            for (var i = 6; i > 0; i--)
            {
                if (currentHour >= forecasts[i])
                {
                    forecast = forecasts[i];
                    break;
                }
            }

            string forecastStr = forecast.ToString();
            if (forecast < 10) { forecastStr = 0 + forecastStr; }

            return forecastStr;
        }
        public static string GetFilename(string cycleStr, string forecast)
        {
            string cycle = int.Parse(cycleStr.Substring(Math.Max(0, cycleStr.Length - 2))).ToString();
            string filename = "gfs.t" + cycle + "z.pgrb2full.0p50.f0" + forecast;
            return filename;
        }

        public static int GetMillibars(double alt)
        {
            if (alt < 710) // 1000mb
            {
                return -9999;
            }
            else if (alt > 54000) // 75mb
            {
                return -9999;
            }

            if (alt >= 710 && alt < 1410)
            {
                return 975;
            }
            else if (alt >= 1410 && alt < 2135)
            {
                return 950;
            }
            else if (alt >= 2135 && alt < 2875)
            {
                return 925;
            }
            else if (alt >= 2875 && alt < 4000)
            {
                return 900;
            }
            else if (alt >= 4000 && alt < 5600)
            {
                return 850;
            }
            else if (alt >= 5600 && alt < 7235)
            {
                return 800;
            }
            else if (alt >= 7235 && alt < 8995)
            {
                return 750;
            }
            else if (alt >= 8995 && alt < 10820)
            {
                return 700;
            }
            else if (alt >= 10820 && alt < 12775)
            {
                return 650;
            }
            else if (alt >= 12775 && alt < 14865)
            {
                return 600;
            }
            else if (alt >= 14865 && alt < 17100)
            {
                return 550;
            }
            else if (alt >= 17100 && alt < 19530)
            {
                return 500;
            }
            else if (alt >= 19530 && alt < 22160)
            {
                return 450;
            }
            else if (alt >= 22160 && alt < 25065)
            {
                return 400;
            }
            else if (alt >= 25065 && alt < 27450)
            {
                return 350;
            }
            else if (alt >= 27450 && alt < 29170)
            {
                return 325;
            }
            else if (alt >= 29170 && alt < 31000)
            {
                return 300;
            }
            else if (alt >= 31000 && alt < 32965)
            {
                return 275;
            }
            else if (alt >= 32965 && alt < 35100)
            {
                return 250;
            }
            else if (alt >= 35100 && alt < 37400)
            {
                return 225;
            }
            else if (alt >= 37400 && alt < 39950)
            {
                return 200;
            }
            else if (alt >= 39950 && alt < 42800)
            {
                return 175;
            }
            else if (alt >= 42800 && alt < 46000)
            {
                return 150;
            }
            else if (alt >= 46000 && alt < 49750)
            {
                return 125;
            }
            else if (alt >= 49750 && alt < 54000)
            {
                return 100;
            }
            return -9999;
        }

        public static int LaunchWGrib2(string loadedFilePath)
        {
            if (Settings.Default.PATHwgrib2 != "")
            {
                PATHwgrib2 = Settings.Default.PATHwgrib2;
            } else {
                Form2 f2 = new Form2();
                f2.ShowDialog(); // Shows Form2
                return 0;
            }

            if (loadedFilePath == "")
            {
                MessageBox.Show("File path cannot be empty", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                return 0;
            }

            ProcessStartInfo startInfo = new ProcessStartInfo
            {
                FileName = "CMD.EXE",
                Arguments = @"/c" + PATHwgrib2 + "\\wgrib2.exe " + loadedFilePath + " -wind_speed " + PATHwgrib2 + "/wind2.grb -wind_dir " + PATHwgrib2 + "/wind2.grb",
                UseShellExecute = false,
                RedirectStandardOutput = true,
                RedirectStandardError = true
            };
            Process p = Process.Start(startInfo);
            p.WaitForExit();
            return 100;
        }

        public static void LaunchWGrib2Again(double lat, double lon)
        {
            if (Settings.Default.PATHwgrib2 != "")
            {
                PATHwgrib2 = Settings.Default.PATHwgrib2;
            }
            else
            {
                Form2 f2 = new Form2();
                f2.ShowDialog(); // Shows Form2
                return;
            }

            if (lon == -9999 || lat == -9999 || alt == -9999)
            {
                MessageBox.Show("No data received from FlightGear", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                return;
            }

            // int levelMB = GetMillibars(alt);

            ProcessStartInfo startInfo = new ProcessStartInfo
            {
                FileName = "CMD.EXE",
                Arguments = @"/c" + PATHwgrib2 + "\\wgrib2.exe " + PATHwgrib2 + "/wind2.grb" + " -s -lon " + lon.ToString() + " " + lat.ToString() + " -match \":(150|200|250|300|400|500|650|850|1000) mb:\"",
                UseShellExecute = false
            };

            Console.WriteLine(startInfo.Arguments);
            startInfo.RedirectStandardOutput = true;
            startInfo.RedirectStandardError = true;
            Process p = Process.Start(startInfo);
            int messageFlag = 0;
            while (!p.StandardOutput.EndOfStream)
            {
                string line = p.StandardOutput.ReadLine();
                Console.WriteLine(line);
                string[] splitLine = line.Split('=');
                switch (messageFlag)
                {
                    case 0:
                        messageFlag = 1;
                        windSpd[0] = double.Parse(splitLine[4].ToString()) * 1.94384;
                        break;
                    case 1:
                        messageFlag = 2;
                        windDir[0] = double.Parse(splitLine[4].ToString());
                        break;
                    case 2:
                        messageFlag = 3;
                        windSpd[1] = double.Parse(splitLine[4].ToString()) * 1.94384;
                        break;
                    case 3:
                        messageFlag = 4;
                        windDir[1] = double.Parse(splitLine[4].ToString());
                        break;
                    case 4:
                        messageFlag = 5;
                        windSpd[2] = double.Parse(splitLine[4].ToString()) * 1.94384;
                        break;
                    case 5:
                        messageFlag = 6;
                        windDir[2] = double.Parse(splitLine[4].ToString());
                        break;
                    case 6:
                        messageFlag = 7;
                        windSpd[3] = double.Parse(splitLine[4].ToString()) * 1.94384;
                        break;
                    case 7:
                        messageFlag = 8;
                        windDir[3] = double.Parse(splitLine[4].ToString());
                        break;
                    case 8:
                        messageFlag = 9;
                        windSpd[4] = double.Parse(splitLine[4].ToString()) * 1.94384;
                        break;
                    case 9:
                        messageFlag = 10;
                        windDir[4] = double.Parse(splitLine[4].ToString());
                        break;
                    case 10:
                        messageFlag = 11;
                        windSpd[5] = double.Parse(splitLine[4].ToString()) * 1.94384;
                        break;
                    case 11:
                        messageFlag = 12;
                        windDir[5] = double.Parse(splitLine[4].ToString());
                        break;
                    case 12:
                        messageFlag = 13;
                        windSpd[6] = double.Parse(splitLine[4].ToString()) * 1.94384;
                        break;
                    case 13:
                        messageFlag = 14;
                        windDir[6] = double.Parse(splitLine[4].ToString());
                        break;
                    case 14:
                        messageFlag = 15;
                        windSpd[7] = double.Parse(splitLine[4].ToString()) * 1.94384;
                        break;
                    case 15:
                        messageFlag = 16;
                        windDir[7] = double.Parse(splitLine[4].ToString());
                        break;
                    case 16:
                        messageFlag = 17;
                        windSpd[8] = double.Parse(splitLine[4].ToString()) * 1.94384;
                        break;
                    case 17:
                        messageFlag = 18;
                        windDir[8] = double.Parse(splitLine[4].ToString());
                        break;
                }
            }
            p.WaitForExit();
            PrepareSendData();
        }

        static void PrepareSendData()
        {
            if (lat != 9999 && lon != 9999)
            {
                SendData(Settings.Default.UDPportIn, lat, lon, windDir[0], windSpd[0], windDir[1], windSpd[1], windDir[2], windSpd[2],
                    windDir[3], windSpd[3], windDir[4], windSpd[4], windDir[5], windSpd[5],
                    windDir[6], windSpd[6], windDir[7], windSpd[7], windDir[8], windSpd[8]);
            }
            else
            {
                MessageBox.Show("Latitude and longitude data not yet received!", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                return;
            }
        }

        static async void OnUdpData(IAsyncResult result)
        {
            UdpClient socket = result.AsyncState as UdpClient;
            IPEndPoint source = new IPEndPoint(0, 0);
            byte[] message = socket.EndReceive(result, ref source);

            string returnData = Encoding.ASCII.GetString(message);
            string[] arrayStr = ParseReceivedData(returnData);
    
            try
            {
                alt = double.Parse(arrayStr[0]);
                double latNew = double.Parse(arrayStr[1]);
                latRounded = (Math.Round(latNew * 2, MidpointRounding.AwayFromZero) / 2);
                double lonNew = double.Parse(arrayStr[2]);
                lonRounded = (Math.Round(lonNew * 2, MidpointRounding.AwayFromZero) / 2);
            }
            catch (Exception ex) {
                string errorString = "Error: " + ex.Message;
                MessageBox.Show(errorString, "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);

                latRounded = -1234;
                lonRounded = -1234;
            }

            _willParseflag = false;

            if (latRounded != -1234 && lonRounded != -1234)
            {
                if ((lat - latRounded) > 0.25 || (lat - latRounded) < -0.25)
                {
                    lat = latRounded;
                    if (_willParseflag == false) { _willParseflag = true; }
                }

                if ((lon - lonRounded) > 0.25 || (lon - lonRounded) < -0.25)
                {
                    lon = lonRounded;
                    if (_willParseflag == false) { _willParseflag = true; }
                }
            }
            else
            {
                MessageBox.Show("Error: invalid latitude / longitude data received", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
            }

            if (_willParseflag == true)
            {
                await Task.Run(() => LaunchWGrib2Again(lat, lon));
            }

            socket.BeginReceive(new AsyncCallback(OnUdpData), socket);
        }

        public static void ReceiveData(int port)
        {
            UdpClient socket = new UdpClient(port);
            socket.BeginReceive(new AsyncCallback(OnUdpData), socket);
            IPEndPoint target = new IPEndPoint(IPAddress.Parse("127.0.0.1"), port);
        }

        static string[] ParseReceivedData(string message)
        {
            string[] messageList = message.Split(',');
            return messageList;
        }

        static void OnUdpDataSend(IAsyncResult result)
        {
            // this is what had been passed into BeginReceive as the second parameter:
            UdpClient socket = result.AsyncState as UdpClient;
            // points towards whoever had sent the message:
            IPEndPoint source = new IPEndPoint(0, 0);
        }


        public static void SendData(int port, double lat, double lon, double d0, double v0, double d1, double v1, double d2, double v2, double d3, double v3, double d4, double v4, double d5, double v5, double d6, double v6, double d7, double v7, double d8, double v8)
        {
            try
            {
                UdpClient socketSend = new UdpClient(port);
                string msg = lat.ToString() + "," + lon.ToString() + "," + d0.ToString() + "," + v0.ToString() + "," + d1.ToString() + "," + v1.ToString() + "," 
                    + d2.ToString() + "," + v2.ToString() + "," + d3.ToString() + "," + v3.ToString() + "," 
                    + d4.ToString() + "," + v4.ToString() + "," + d5.ToString() + "," + v5.ToString() + "," 
                    + d6.ToString() + "," + v6.ToString() + "," + d7.ToString() + "," + v7.ToString() + "," 
                    + d8.ToString() + "," + v8.ToString();

                msg += "\n";
                byte[] payload = Encoding.ASCII.GetBytes(msg);

                IPEndPoint target = new IPEndPoint(IPAddress.Parse("127.0.0.1"), port);
                socketSend.Connect(target); ;

                //Console.WriteLine("Sending " + payload.Length + " bytes to " + target);
                string returnData2 = Encoding.ASCII.GetString(payload);
                //Console.WriteLine(returnData2);

                socketSend.Send(payload, payload.Length);
                socketSend.Close();
            }
            catch
            {
                MessageBox.Show("Unable to start connection - socket already in use!", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
            }
            }
    }
}
