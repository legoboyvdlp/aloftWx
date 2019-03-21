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
            this.GetCycleButton = new System.Windows.Forms.Button();
            this.progressBar1 = new System.Windows.Forms.ProgressBar();
            this.label1 = new System.Windows.Forms.Label();
            this.curProgress = new System.Windows.Forms.Label();
            this.urlBox = new System.Windows.Forms.TextBox();
            this.SuspendLayout();
            // 
            // GetCycleButton
            // 
            this.GetCycleButton.Location = new System.Drawing.Point(11, 25);
            this.GetCycleButton.Margin = new System.Windows.Forms.Padding(2);
            this.GetCycleButton.Name = "GetCycleButton";
            this.GetCycleButton.Size = new System.Drawing.Size(97, 28);
            this.GetCycleButton.TabIndex = 2;
            this.GetCycleButton.Text = "Fetch Cycle";
            this.GetCycleButton.UseVisualStyleBackColor = true;
            this.GetCycleButton.Click += new System.EventHandler(this.Button1_ClickAsync);
            // 
            // progressBar1
            // 
            this.progressBar1.Location = new System.Drawing.Point(12, 73);
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
            // AloftWxForm
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(435, 177);
            this.Controls.Add(this.urlBox);
            this.Controls.Add(this.curProgress);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.progressBar1);
            this.Controls.Add(this.GetCycleButton);
            this.Margin = new System.Windows.Forms.Padding(2);
            this.Name = "AloftWxForm";
            this.Text = "AloftWx";
            this.Load += new System.EventHandler(this.Form1_Load);
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion
        private System.Windows.Forms.Button GetCycleButton;
        private System.Windows.Forms.ProgressBar progressBar1;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.Label curProgress;
        private System.Windows.Forms.TextBox urlBox;
    }
}

