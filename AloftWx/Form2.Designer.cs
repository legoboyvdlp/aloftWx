namespace AloftWx
{
    partial class Form2
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
            this.label1 = new System.Windows.Forms.Label();
            this.wgrib2pathbox = new System.Windows.Forms.TextBox();
            this.button1 = new System.Windows.Forms.Button();
            this.button2 = new System.Windows.Forms.Button();
            this.udpButton = new System.Windows.Forms.Button();
            this.udpBox = new System.Windows.Forms.TextBox();
            this.label2 = new System.Windows.Forms.Label();
            this.udpButtonIn = new System.Windows.Forms.Button();
            this.udpBoxIn = new System.Windows.Forms.TextBox();
            this.label3 = new System.Windows.Forms.Label();
            this.SuspendLayout();
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(14, 16);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(72, 13);
            this.label1.TabIndex = 0;
            this.label1.Text = "Path to wgrib:";
            // 
            // wgrib2pathbox
            // 
            this.wgrib2pathbox.Location = new System.Drawing.Point(92, 12);
            this.wgrib2pathbox.Name = "wgrib2pathbox";
            this.wgrib2pathbox.Size = new System.Drawing.Size(212, 20);
            this.wgrib2pathbox.TabIndex = 1;
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(310, 12);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(22, 22);
            this.button1.TabIndex = 2;
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.Button1_Click);
            // 
            // button2
            // 
            this.button2.Location = new System.Drawing.Point(29, 161);
            this.button2.Name = "button2";
            this.button2.Size = new System.Drawing.Size(274, 22);
            this.button2.TabIndex = 3;
            this.button2.Text = "Save";
            this.button2.UseVisualStyleBackColor = true;
            this.button2.Click += new System.EventHandler(this.button2_Click);
            // 
            // udpButton
            // 
            this.udpButton.Location = new System.Drawing.Point(310, 41);
            this.udpButton.Name = "udpButton";
            this.udpButton.Size = new System.Drawing.Size(22, 22);
            this.udpButton.TabIndex = 6;
            this.udpButton.UseVisualStyleBackColor = true;
            this.udpButton.Click += new System.EventHandler(this.udpButton_Click);
            // 
            // udpBox
            // 
            this.udpBox.Location = new System.Drawing.Point(92, 41);
            this.udpBox.Name = "udpBox";
            this.udpBox.Size = new System.Drawing.Size(212, 20);
            this.udpBox.TabIndex = 5;
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(14, 44);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(78, 13);
            this.label2.TabIndex = 4;
            this.label2.Text = "UDP port (out):";
            // 
            // udpButtonIn
            // 
            this.udpButtonIn.Location = new System.Drawing.Point(310, 71);
            this.udpButtonIn.Name = "udpButtonIn";
            this.udpButtonIn.Size = new System.Drawing.Size(22, 22);
            this.udpButtonIn.TabIndex = 9;
            this.udpButtonIn.UseVisualStyleBackColor = true;
            this.udpButtonIn.Click += new System.EventHandler(this.udpButtonIn_Click);
            // 
            // udpBoxIn
            // 
            this.udpBoxIn.Location = new System.Drawing.Point(92, 71);
            this.udpBoxIn.Name = "udpBoxIn";
            this.udpBoxIn.Size = new System.Drawing.Size(212, 20);
            this.udpBoxIn.TabIndex = 8;
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(14, 74);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(71, 13);
            this.label3.TabIndex = 7;
            this.label3.Text = "UDP port (in):";
            // 
            // Form2
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(334, 186);
            this.Controls.Add(this.udpButtonIn);
            this.Controls.Add(this.udpBoxIn);
            this.Controls.Add(this.label3);
            this.Controls.Add(this.udpButton);
            this.Controls.Add(this.udpBox);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.button2);
            this.Controls.Add(this.button1);
            this.Controls.Add(this.wgrib2pathbox);
            this.Controls.Add(this.label1);
            this.Name = "Form2";
            this.Text = "Options";
            this.Load += new System.EventHandler(this.Form2_Load);
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.TextBox wgrib2pathbox;
        private System.Windows.Forms.Button button1;
        private System.Windows.Forms.Button button2;
        private System.Windows.Forms.Button udpButton;
        private System.Windows.Forms.TextBox udpBox;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.Button udpButtonIn;
        private System.Windows.Forms.TextBox udpBoxIn;
        private System.Windows.Forms.Label label3;
    }
}