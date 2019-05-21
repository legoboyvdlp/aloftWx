using System.Drawing;

namespace AloftWx
{
    partial class AloftWxForm
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.components = new System.ComponentModel.Container();
            this.GetCycleButton = new System.Windows.Forms.Button();
            this.progressBar1 = new System.Windows.Forms.ProgressBar();
            this.label1 = new System.Windows.Forms.Label();
            this.curProgress = new System.Windows.Forms.Label();
            this.urlBox = new System.Windows.Forms.TextBox();
            this.menuStrip1 = new System.Windows.Forms.MenuStrip();
            this.fileToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.quitToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.optionsToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.connectToFGToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.button1 = new System.Windows.Forms.Button();
            this.loadedfile = new System.Windows.Forms.Label();
            this.loadedfileBox = new System.Windows.Forms.TextBox();
            this.button2 = new System.Windows.Forms.Button();
            this.button3 = new System.Windows.Forms.Button();
            this.progressBar2 = new System.Windows.Forms.ProgressBar();
            this.contextMenuStrip1 = new System.Windows.Forms.ContextMenuStrip(this.components);
            this.label2 = new System.Windows.Forms.Label();
            this.statusLabel = new System.Windows.Forms.Label();
            this.menuStrip1.SuspendLayout();
            this.SuspendLayout();
            // 
            // GetCycleButton
            // 
            this.GetCycleButton.Location = new System.Drawing.Point(11, 25);
            this.GetCycleButton.Margin = new System.Windows.Forms.Padding(2);
            this.GetCycleButton.Name = "GetCycleButton";
            this.GetCycleButton.Size = new System.Drawing.Size(97, 28);
            this.GetCycleButton.TabIndex = 2;
            this.GetCycleButton.Text = "Download Cycle";
            this.GetCycleButton.UseVisualStyleBackColor = true;
            this.GetCycleButton.Click += new System.EventHandler(this.Button1_ClickAsync);
            // 
            // progressBar1
            // 
            this.progressBar1.Location = new System.Drawing.Point(12, 58);
            this.progressBar1.Name = "progressBar1";
            this.progressBar1.Size = new System.Drawing.Size(411, 28);
            this.progressBar1.TabIndex = 3;
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Font = new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.label1.Location = new System.Drawing.Point(145, 33);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(36, 13);
            this.label1.TabIndex = 8;
            this.label1.Text = "URL:";
            // 
            // curProgress
            // 
            this.curProgress.AutoSize = true;
            this.curProgress.Location = new System.Drawing.Point(297, 155);
            this.curProgress.Name = "curProgress";
            this.curProgress.Size = new System.Drawing.Size(0, 13);
            this.curProgress.TabIndex = 13;
            // 
            // urlBox
            // 
            this.urlBox.Location = new System.Drawing.Point(197, 30);
            this.urlBox.Name = "urlBox";
            this.urlBox.Size = new System.Drawing.Size(226, 20);
            this.urlBox.TabIndex = 14;
            // 
            // menuStrip1
            // 
            this.menuStrip1.BackColor = System.Drawing.SystemColors.ControlLight;
            this.menuStrip1.Items.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.fileToolStripMenuItem,
            this.optionsToolStripMenuItem,
            this.connectToFGToolStripMenuItem});
            this.menuStrip1.Location = new System.Drawing.Point(0, 0);
            this.menuStrip1.Name = "menuStrip1";
            this.menuStrip1.Size = new System.Drawing.Size(435, 24);
            this.menuStrip1.TabIndex = 16;
            this.menuStrip1.Text = "menuStrip1";
            // 
            // fileToolStripMenuItem
            // 
            this.fileToolStripMenuItem.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.quitToolStripMenuItem});
            this.fileToolStripMenuItem.Name = "fileToolStripMenuItem";
            this.fileToolStripMenuItem.Size = new System.Drawing.Size(37, 20);
            this.fileToolStripMenuItem.Text = "File";
            // 
            // quitToolStripMenuItem
            // 
            this.quitToolStripMenuItem.Name = "quitToolStripMenuItem";
            this.quitToolStripMenuItem.ShortcutKeys = ((System.Windows.Forms.Keys)(((System.Windows.Forms.Keys.Control | System.Windows.Forms.Keys.Shift) 
            | System.Windows.Forms.Keys.Q)));
            this.quitToolStripMenuItem.Size = new System.Drawing.Size(172, 22);
            this.quitToolStripMenuItem.Text = "Quit";
            this.quitToolStripMenuItem.Click += new System.EventHandler(this.QuitToolStripMenuItem_Click);
            // 
            // optionsToolStripMenuItem
            // 
            this.optionsToolStripMenuItem.Name = "optionsToolStripMenuItem";
            this.optionsToolStripMenuItem.Size = new System.Drawing.Size(61, 20);
            this.optionsToolStripMenuItem.Text = "Options";
            this.optionsToolStripMenuItem.Click += new System.EventHandler(this.OptionsToolStripMenuItem_Click_1);
            // 
            // connectToFGToolStripMenuItem
            // 
            this.connectToFGToolStripMenuItem.Name = "connectToFGToolStripMenuItem";
            this.connectToFGToolStripMenuItem.Size = new System.Drawing.Size(113, 20);
            this.connectToFGToolStripMenuItem.Text = "Open Connection";
            this.connectToFGToolStripMenuItem.Click += new System.EventHandler(this.ConnectToFGToolStripMenuItem_Click);
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(12, 138);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(96, 46);
            this.button1.TabIndex = 17;
            this.button1.Text = "Parse loaded .grib2 file";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.Button1_Click);
            // 
            // loadedfile
            // 
            this.loadedfile.AutoSize = true;
            this.loadedfile.Font = new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.loadedfile.Location = new System.Drawing.Point(113, 162);
            this.loadedfile.Name = "loadedfile";
            this.loadedfile.Size = new System.Drawing.Size(74, 13);
            this.loadedfile.TabIndex = 18;
            this.loadedfile.Text = "Loaded file:";
            // 
            // loadedfileBox
            // 
            this.loadedfileBox.Location = new System.Drawing.Point(193, 159);
            this.loadedfileBox.Name = "loadedfileBox";
            this.loadedfileBox.Size = new System.Drawing.Size(226, 20);
            this.loadedfileBox.TabIndex = 19;
            // 
            // button2
            // 
            this.button2.Location = new System.Drawing.Point(12, 93);
            this.button2.Name = "button2";
            this.button2.Size = new System.Drawing.Size(95, 41);
            this.button2.TabIndex = 20;
            this.button2.Text = "Load downloaded file";
            this.button2.UseVisualStyleBackColor = true;
            this.button2.Click += new System.EventHandler(this.Button2_Click);
            // 
            // button3
            // 
            this.button3.Location = new System.Drawing.Point(116, 93);
            this.button3.Name = "button3";
            this.button3.Size = new System.Drawing.Size(95, 41);
            this.button3.TabIndex = 21;
            this.button3.Text = "Load other grib2 file";
            this.button3.UseVisualStyleBackColor = true;
            this.button3.Click += new System.EventHandler(this.Button3_Click);
            // 
            // progressBar2
            // 
            this.progressBar2.Location = new System.Drawing.Point(12, 190);
            this.progressBar2.Name = "progressBar2";
            this.progressBar2.Size = new System.Drawing.Size(411, 28);
            this.progressBar2.TabIndex = 22;
            // 
            // contextMenuStrip1
            // 
            this.contextMenuStrip1.Name = "contextMenuStrip1";
            this.contextMenuStrip1.Size = new System.Drawing.Size(61, 4);
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Font = new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.label2.Location = new System.Drawing.Point(12, 227);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(47, 13);
            this.label2.TabIndex = 24;
            this.label2.Text = "Status:";
            // 
            // statusLabel
            // 
            this.statusLabel.AutoSize = true;
            this.statusLabel.Font = new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.statusLabel.Location = new System.Drawing.Point(65, 227);
            this.statusLabel.Name = "statusLabel";
            this.statusLabel.Size = new System.Drawing.Size(67, 13);
            this.statusLabel.TabIndex = 25;
            this.statusLabel.Text = "Not Running";
            // 
            // AloftWxForm
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(435, 249);
            this.Controls.Add(this.statusLabel);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.progressBar2);
            this.Controls.Add(this.button3);
            this.Controls.Add(this.button2);
            this.Controls.Add(this.loadedfileBox);
            this.Controls.Add(this.loadedfile);
            this.Controls.Add(this.button1);
            this.Controls.Add(this.urlBox);
            this.Controls.Add(this.curProgress);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.progressBar1);
            this.Controls.Add(this.GetCycleButton);
            this.Controls.Add(this.menuStrip1);
            this.Margin = new System.Windows.Forms.Padding(2);
            this.Name = "AloftWxForm";
            this.Text = "AloftWx";
            this.menuStrip1.ResumeLayout(false);
            this.menuStrip1.PerformLayout();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion
        private System.Windows.Forms.Button GetCycleButton;
        private System.Windows.Forms.ProgressBar progressBar1;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.Label curProgress;
        private System.Windows.Forms.TextBox urlBox;
        private System.Windows.Forms.MenuStrip menuStrip1;
        private System.Windows.Forms.ToolStripMenuItem fileToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem quitToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem optionsToolStripMenuItem;
        private System.Windows.Forms.Button button1;
        private System.Windows.Forms.Label loadedfile;
        private System.Windows.Forms.TextBox loadedfileBox;
        private System.Windows.Forms.Button button2;
        private System.Windows.Forms.Button button3;
        private System.Windows.Forms.ToolStripMenuItem connectToFGToolStripMenuItem;
        private System.Windows.Forms.ProgressBar progressBar2;
        private System.Windows.Forms.ContextMenuStrip contextMenuStrip1;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.Label statusLabel;
    }
}