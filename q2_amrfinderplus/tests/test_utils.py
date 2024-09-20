from unittest.mock import call, patch

from qiime2.plugin.testing import TestPluginBase

from q2_amrfinderplus.utils import EXTERNAL_CMD_WARNING, run_command


class TestRunCommand(TestPluginBase):
    package = "q2_amrfinderplus.tests"

    @patch("subprocess.run")
    @patch("builtins.print")
    def test_run_command_verbose(self, mock_print, mock_subprocess_run):
        # Mock command and working directory
        cmd = ["echo", "Hello"]
        cwd = "/test/directory"

        # Run the function with verbose=True
        run_command(cmd, cwd=cwd, verbose=True)

        # Check if subprocess.run was called with the correct arguments
        mock_subprocess_run.assert_called_once_with(cmd, check=True, cwd=cwd)

        # Check if the correct print statements were called
        mock_print.assert_has_calls(
            [
                call(EXTERNAL_CMD_WARNING),
                call("\nCommand:", end=" "),
                call("echo Hello", end="\n\n"),
            ]
        )

    @patch("subprocess.run")
    @patch("builtins.print")
    def test_run_command_non_verbose(self, mock_print, mock_subprocess_run):
        # Mock command and working directory
        cmd = ["echo", "Hello"]
        cwd = "/test/directory"

        # Run the function with verbose=False
        run_command(cmd, cwd=cwd, verbose=False)

        # Check if subprocess.run was called with the correct arguments
        mock_subprocess_run.assert_called_once_with(cmd, check=True, cwd=cwd)

        # Ensure no print statements were made
        mock_print.assert_not_called()
