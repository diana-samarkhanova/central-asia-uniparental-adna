from __future__ import annotations

import contextlib
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from urllib.error import HTTPError, URLError

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import fetch_sources as fetch


class Response(io.BytesIO):
    def geturl(self):
        return "https://example.org/frozen"


class FetchSourcesTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.body = b"the exact frozen metadata\n"
        self.source = {"key": "aadr", "release": "test frozen version",
                       "expected_path": "data/raw/aadr/input.anno", "download_url": "https://example.org/frozen",
                       "sha256": hashlib.sha256(self.body).hexdigest()}
        self.destination = self.root / self.source["expected_path"]

    def acquire(self, **kwargs):
        return fetch.acquire(self.source, self.root, **kwargs)

    def assert_no_parts(self):
        self.assertEqual(list(self.root.rglob("*.part")), [])

    def test_download_verifies_hash_and_installs_once(self):
        with patch.object(fetch, "urlopen", return_value=Response(self.body)) as request:
            result = self.acquire()
        self.assertEqual(result["status"], "READY")
        self.assertEqual(self.destination.read_bytes(), self.body)
        self.assertEqual(result["bytes"], len(self.body))
        self.assertEqual(request.call_args.kwargs["timeout"], 30)
        with patch.object(fetch, "urlopen") as request:
            self.assertEqual(self.acquire()["status"], "READY")
            request.assert_not_called()
        self.assert_no_parts()

    def test_mismatching_download_is_not_installed_or_retried(self):
        with patch.object(fetch, "urlopen", return_value=Response(b"wrong release")) as request:
            result = self.acquire()
        self.assertEqual(result["status"], "HASH_MISMATCH")
        self.assertEqual(request.call_count, 1)
        self.assertFalse(self.destination.exists())
        self.assert_no_parts()

    def test_existing_conflict_is_unchanged(self):
        self.destination.parent.mkdir(parents=True)
        self.destination.write_bytes(b"keep these bytes")
        with patch.object(fetch, "urlopen") as request:
            self.assertEqual(self.acquire()["status"], "CONFLICT")
            request.assert_not_called()
        self.assertEqual(self.destination.read_bytes(), b"keep these bytes")

    def test_local_import_and_missing_or_mismatching_local(self):
        retained = self.root / "retained.csv"
        self.assertEqual(self.acquire(local=retained)["status"], "ERROR")
        retained.write_bytes(b"newer version")
        self.assertEqual(self.acquire(local=retained)["status"], "HASH_MISMATCH")
        retained.write_bytes(self.body)
        with patch.object(fetch, "urlopen") as request:
            self.assertEqual(self.acquire(local=retained)["status"], "READY")
            request.assert_not_called()
        self.assertEqual(retained.read_bytes(), self.body)
        self.assert_no_parts()

    def test_blocked_amtdb_does_not_use_live_export(self):
        self.source.update(key="amtdb", download_url="")
        with patch.object(fetch, "urlopen") as request:
            result = self.acquire()
            request.assert_not_called()
        self.assertEqual(result["status"], "BLOCKED")
        self.assertIn("--local amtdb=", result["message"])
        self.assertIn("cannot substitute", result["message"])

    def test_offline_absent_source(self):
        with patch.object(fetch, "urlopen") as request:
            self.assertEqual(self.acquire(offline=True)["status"], "MISSING")
            request.assert_not_called()

    def test_transient_network_errors_retry_then_succeed(self):
        for error in (HTTPError("u", 503, "Unavailable", {}, None), URLError(TimeoutError("timed out"))):
            with self.subTest(error=type(error).__name__):
                with patch.object(fetch, "urlopen", side_effect=[error, Response(self.body)]) as request, patch.object(fetch.time, "sleep"):
                    result = self.acquire()
                self.assertEqual(result["status"], "READY")
                self.assertEqual(result["attempts"], 2)
                self.assertEqual(request.call_count, 2)
                self.destination.unlink()

    def test_permanent_failure_never_retries(self):
        for error in (HTTPError("u", 404, "Not found", {}, None), URLError("access restricted")):
            with patch.object(fetch, "urlopen", side_effect=error) as request:
                self.assertEqual(self.acquire()["status"], "ERROR")
            self.assertEqual(request.call_count, 1)
        self.assert_no_parts()

    def test_retry_budget_is_bounded(self):
        with patch.object(fetch, "urlopen", side_effect=TimeoutError("timeout")) as request, patch.object(fetch.time, "sleep"):
            result = self.acquire(retries=2)
        self.assertEqual(result["status"], "ERROR")
        self.assertEqual(request.call_count, 3)

    def test_interruption_removes_partial_download(self):
        class Interrupted(Response):
            def read(self, size):
                if self.tell():
                    raise KeyboardInterrupt
                return super().read(4)
        with patch.object(fetch, "urlopen", return_value=Interrupted(self.body)):
            result = self.acquire()
        self.assertEqual(result["status"], "INTERRUPTED")
        self.assertFalse(self.destination.exists())
        self.assert_no_parts()

    def test_midstream_retry_discards_partial_bytes(self):
        class Interrupted(Response):
            def read(self, size):
                if self.tell():
                    raise ConnectionResetError("connection reset")
                return super().read(4)
        with patch.object(fetch, "urlopen", side_effect=[Interrupted(self.body), Response(self.body)]), patch.object(fetch.time, "sleep"):
            result = self.acquire()
        self.assertEqual(result["status"], "READY")
        self.assertEqual(self.destination.read_bytes(), self.body)
        self.assert_no_parts()

    def test_concurrent_conflicting_destination_is_never_overwritten(self):
        real_link = fetch.os.link
        def race(source, destination):
            destination.write_bytes(b"other process")
            return real_link(source, destination)
        with patch.object(fetch, "urlopen", return_value=Response(self.body)), patch.object(fetch.os, "link", side_effect=race):
            result = self.acquire()
        self.assertEqual(result["status"], "CONFLICT")
        self.assertEqual(self.destination.read_bytes(), b"other process")
        self.assert_no_parts()

    def test_redirect_downgrade_and_oversized_file_are_rejected(self):
        response = Response(self.body)
        response.geturl = lambda: "http://example.org/insecure"
        with patch.object(fetch, "urlopen", return_value=response):
            self.assertEqual(self.acquire()["status"], "ERROR")
        with patch.object(fetch, "urlopen", return_value=Response(self.body)), patch.object(fetch, "MAX_BYTES", 4):
            self.assertEqual(self.acquire()["status"], "ERROR")
        self.assertFalse(self.destination.exists())
        self.assert_no_parts()

    def test_registry_has_only_known_frozen_hashes_and_immutable_urls(self):
        sources = fetch.load_sources(ROOT / "data/SOURCES.tsv")
        self.assertEqual(set(sources), set(fetch.EXPECTED))
        self.assertEqual(sources["amtdb"]["download_url"], "")
        self.assertIn("/13994518", sources["aadr"]["download_url"])
        self.assertIn("bc770a59ace8cd4c042c6f903d620d93ee751eb0", sources["aychr"]["download_url"])

    def test_cli_selection_report_and_incomplete_exit(self):
        report_path = self.root / "report.json"
        with contextlib.redirect_stdout(io.StringIO()) as output:
            code = fetch.main(["--resource", "amtdb", "--output-root", str(self.root), "--report", str(report_path)])
        report = json.loads(output.getvalue())
        self.assertEqual(code, 1)
        self.assertEqual(report["selected_resources"], ["amtdb"])
        self.assertFalse(report["all_three_sources_ready"])
        self.assertEqual(report["sources"][0]["status"], "BLOCKED")
        self.assertEqual(report, json.loads(report_path.read_text()))

    def test_cli_report_cannot_overwrite_input(self):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as caught:
            fetch.main(["--report", str(ROOT / "data/SOURCES.tsv")])
        self.assertEqual(caught.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
