from unittest.mock import patch

from qiime2.plugin.testing import TestPluginBase

from q2_amrfinderplus.database import fetch_amrfinderplus_db, run_amrfinder_fetch


class TestFetchAMRFinderPlusDB(TestPluginBase):
    package = "q2_amrfinderplus.tests"

    @patch("q2_amrfinderplus.database.run_amrfinder_fetch")
    @patch("q2_amrfinderplus.database._copy_all")
    def test_fetch_amrfinderplus_db(self, mock_run_amrfinder_u, mock__copy_all):
        fetch_amrfinderplus_db()

    @patch("q2_amrfinderplus.database.run_command")
    def test_run_amrfinder_u(self, mock_run_command):
        run_amrfinder_fetch()
        mock_run_command.assert_called_once_with(
            ["amrfinder", "-u"],
            verbose=True,
        )
