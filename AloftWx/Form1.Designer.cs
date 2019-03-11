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
            this.GetCycleButton.Click += new System.EventHandler(this.Button1_Click);
            // 
            // AloftWxForm
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(496, 473);
            this.Controls.Add(this.GetCycleButton);
            this.Margin = new System.Windows.Forms.Padding(2);
            this.Name = "AloftWxForm";
            this.Text = "AloftWx";
            this.Load += new System.EventHandler(this.Form1_Load);
            this.ResumeLayout(false);

        }

        #endregion
        private System.Windows.Forms.Button GetCycleButton;
    }
}

