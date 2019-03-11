using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace AloftWx
{
    static class Program
    {
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
            for (var i = 3; i > 0; i--)
            {
                if (currentHour >= cycles[i])
                {
                    curCycle = cycles[i];
                    break;
                }
            }

            // if the cycle is 0, add an extra 0
            string curCycleStr = curCycle.ToString();
            curCycleStr = 0 + curCycleStr;

            // if the month is less than 10, add an extra 0
            string currentMonthStr = currentMonth.ToString();
            if (currentMonth < 10) { currentMonthStr = 0 + currentMonth.ToString(); }

            return currentYear.ToString() + currentMonthStr + currentDay.ToString() + curCycleStr;
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
            string filename = "WAFS_blended_"  + cycleStr + "f" + forecast + ".grib2";
            return filename;
        }
    }
}
