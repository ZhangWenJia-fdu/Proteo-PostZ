using System;
using System.Diagnostics;
using System.IO;
using System.Net.Sockets;
using System.Text;
using System.Threading;
using System.Windows.Forms;

class Launcher {
  const string AppTitle = "ProteoPostZ Formal V2.0";
  const string Url = "http://127.0.0.1:3840/";
  const string Host = "127.0.0.1";
  const int Port = 3840;

  static bool IsPortOpen() {
    try {
      using (var client = new TcpClient()) {
        var result = client.BeginConnect(Host, Port, null, null);
        bool ok = result.AsyncWaitHandle.WaitOne(TimeSpan.FromMilliseconds(700));
        if (!ok) return false;
        client.EndConnect(result);
        return true;
      }
    } catch { return false; }
  }

  static void OpenBrowser() {
    try {
      Process.Start("explorer.exe", Url);
      return;
    } catch { }
    try {
      var psi = new ProcessStartInfo();
      psi.FileName = "cmd.exe";
      psi.Arguments = "/c start \"\" \"" + Url + "\"";
      psi.CreateNoWindow = true;
      psi.UseShellExecute = false;
      Process.Start(psi);
      return;
    } catch { }
    try {
      Process.Start("rundll32.exe", "url.dll,FileProtocolHandler " + Url);
    } catch { }
  }

  static string ExtractLineValue(string text, string prefix, bool addSpaces) {
    using (var reader = new StringReader(text)) {
      string line;
      while ((line = reader.ReadLine()) != null) {
        if (line.StartsWith(prefix)) {
          string value = line.Substring(prefix.Length).Trim();
          return addSpaces ? value.Replace(",", ", ") : value;
        }
      }
    }
    return "";
  }

  static int RunRScript(string rscript, string root, string rlibs, string script, string args, string launcherLog, out string output) {
    output = "";
    var psi = new ProcessStartInfo();
    psi.FileName = rscript;
    psi.Arguments = "\"" + script + "\"" + args;
    psi.WorkingDirectory = root;
    psi.UseShellExecute = false;
    psi.CreateNoWindow = true;
    psi.RedirectStandardOutput = true;
    psi.RedirectStandardError = true;
    psi.EnvironmentVariables["R_LIBS_USER"] = rlibs;

    var proc = new Process();
    var buffer = new StringBuilder();
    proc.StartInfo = psi;
    proc.OutputDataReceived += (s, e) => { if (e.Data != null) buffer.AppendLine(e.Data); };
    proc.ErrorDataReceived += (s, e) => { if (e.Data != null) buffer.AppendLine(e.Data); };
    proc.Start();
    proc.BeginOutputReadLine();
    proc.BeginErrorReadLine();
    bool exited = proc.WaitForExit(1200000);
    if (!exited) {
      try { proc.Kill(); } catch { }
      try { File.AppendAllText(launcherLog, DateTime.Now.ToString("yyyy-MM-dd HH:mm:ss") + " R helper timed out." + Environment.NewLine); } catch { }
      output = "R helper timed out.";
      return 99;
    }
    proc.WaitForExit();
    output = buffer.ToString();
    try {
      File.AppendAllText(launcherLog,
        DateTime.Now.ToString("yyyy-MM-dd HH:mm:ss") + Environment.NewLine +
        "Command: " + psi.FileName + " " + psi.Arguments + Environment.NewLine +
        "ExitCode: " + proc.ExitCode + Environment.NewLine +
        output + Environment.NewLine);
    } catch { }
    return proc.ExitCode;
  }

  static bool RunDependencyCheck(string rscript, string root, string rlibs) {
    string checkScript = Path.Combine(root, "check_dependencies.R");
    string logDir = Path.Combine(root, "logs");
    string launcherLog = Path.Combine(logDir, "dependency_check_launcher.log");

    if (!File.Exists(checkScript)) {
      MessageBox.Show("Cannot find check_dependencies.R next to the launcher.\n\nThe app was not started.", AppTitle);
      return false;
    }

    try { Directory.CreateDirectory(logDir); } catch { }

    try {
      string output;
      int exitCode = RunRScript(rscript, root, rlibs, checkScript, "", launcherLog, out output);
      if (exitCode == 0) return true;

      string missing = ExtractLineValue(output, "MISSING_PACKAGES:", true);
      string missingCran = ExtractLineValue(output, "MISSING_CRAN:", true);
      string missingBioc = ExtractLineValue(output, "MISSING_BIOC:", true);
      string missingSystem = ExtractLineValue(output, "MISSING_SYSTEM:", true);

      if (missing.Length > 0) {
        string message = "The app cannot start because required R packages are missing:\n\n" + missing;
        if (missingCran.Length > 0) message += "\n\nCRAN packages:\n" + missingCran;
        if (missingBioc.Length > 0) message += "\n\nBioconductor packages:\n" + missingBioc;
        if (missingSystem.Length > 0) message += "\n\nThese system/base dependencies are not suitable for automatic installation and must be fixed manually:\n" + missingSystem;
        message += "\n\nDo you want to try installing the missing CRAN/Bioconductor R packages into this app's portable R library now?\n\nNo global user R library will be used intentionally. Network access is required.";

        DialogResult choice = MessageBox.Show(message, AppTitle, MessageBoxButtons.YesNo, MessageBoxIcon.Warning);
        if (choice != DialogResult.Yes) {
          MessageBox.Show("Startup stopped because dependency installation was canceled. No packages were installed.", AppTitle);
          return false;
        }

        string installOutput;
        int installExitCode = RunRScript(rscript, root, rlibs, checkScript, " --install-missing", launcherLog, out installOutput);
        if (installExitCode == 0) return true;

        string installError = ExtractLineValue(installOutput, "INSTALL_FAILED:", true);
        if (installError.Length == 0) installError = installOutput.Trim();
        MessageBox.Show(
          "The missing R packages could not be installed. The app was not started.\n\n" +
          installError +
          "\n\nPlease check your network connection and logs\\dependency_install.log.",
          AppTitle);
        return false;
      }

      MessageBox.Show(
        "The app cannot start because the startup dependency check failed.\n\nDetails were written to logs\\dependency_check.log and logs\\dependency_check_launcher.log.",
        AppTitle);
      return false;
    } catch (Exception ex) {
      try { File.AppendAllText(launcherLog, DateTime.Now.ToString("yyyy-MM-dd HH:mm:ss") + " Dependency check failed: " + ex.Message + Environment.NewLine); } catch { }
      MessageBox.Show(
        "The app cannot start because the startup dependency check could not run:\n" +
        ex.Message,
        AppTitle);
      return false;
    }
  }

  [STAThread]
  static void Main() {
    string root = AppDomain.CurrentDomain.BaseDirectory;
    string portableR = Path.Combine(root, "portable", "R-4.5.1", "bin", "x64", "Rscript.exe");
    string rscript = File.Exists(portableR) ? portableR : "Rscript.exe";
    string rlibs = Path.Combine(root, "portable", "Rlibs");
    string appScript = Path.Combine(root, "run_app.R");
    string stdoutLog = Path.Combine(root, "shiny_stdout.log");
    string stderrLog = Path.Combine(root, "shiny_stderr.log");

    if (!File.Exists(appScript)) {
      MessageBox.Show("Cannot find run_app.R next to the launcher.", AppTitle);
      return;
    }

    try {
      if (!IsPortOpen()) {
        if (!RunDependencyCheck(rscript, root, rlibs)) return;

        var psi = new ProcessStartInfo();
        psi.FileName = rscript;
        psi.Arguments = "\"" + appScript + "\"";
        psi.WorkingDirectory = root;
        psi.UseShellExecute = false;
        psi.CreateNoWindow = true;
        psi.RedirectStandardOutput = true;
        psi.RedirectStandardError = true;
        psi.EnvironmentVariables["R_LIBS_USER"] = rlibs;
        var proc = new Process();
        proc.StartInfo = psi;
        proc.OutputDataReceived += (s, e) => { if (e.Data != null) File.AppendAllText(stdoutLog, e.Data + Environment.NewLine); };
        proc.ErrorDataReceived += (s, e) => { if (e.Data != null) File.AppendAllText(stderrLog, e.Data + Environment.NewLine); };
        proc.Start();
        proc.BeginOutputReadLine();
        proc.BeginErrorReadLine();
      }

      for (int i = 0; i < 20; i++) {
        if (IsPortOpen()) break;
        Thread.Sleep(500);
      }

      OpenBrowser();
    } catch (Exception ex) {
      MessageBox.Show("Failed to start the app:\n" + ex.Message + "\n\nPlease open " + Url + " manually if the app is already running.", AppTitle);
    }
  }
}
