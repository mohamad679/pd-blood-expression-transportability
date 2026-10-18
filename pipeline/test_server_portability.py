"""Portability checks using synthetic /proc fixtures; no study records read."""
import unittest
from pathlib import Path
from unittest.mock import patch
from . import verify_permutation_executor as verifier
from . import permutation_executor as executor


class LinuxCounters(unittest.TestCase):
    def test_linux_memory_uses_kibibytes(self):
        with patch.object(verifier.platform, 'system', return_value='Linux'), \
             patch.object(Path, 'read_text', return_value='MemTotal: 33554432 kB\nSwapTotal: 0 kB\n'):
            self.assertEqual(verifier.physical_memory(), 32 * 1024**3)

    def test_linux_swap_distinct_from_page_writeout(self):
        with patch.object(verifier.platform, 'system', return_value='Linux'), \
             patch.object(Path, 'read_text', return_value='pswpout 19\npgpgout 423000\npswpin 3\n'):
            self.assertEqual(verifier.swap_counters(), {'pswpout': 19, 'pgpgout': 423000})

    def test_missing_linux_counter_is_not_zero(self):
        with patch.object(verifier.platform, 'system', return_value='Linux'), \
             patch.object(Path, 'read_text', return_value='pgpgout 0\n'):
            self.assertEqual(verifier.swap_counters(), {})

    def test_missing_memory_cannot_pass(self):
        with patch.object(verifier.platform, 'system', return_value='Linux'), \
             patch.object(Path, 'read_text', side_effect=FileNotFoundError):
            self.assertIsNone(verifier.physical_memory())

    def test_mac_counters_retained(self):
        with patch.object(verifier.platform, 'system', return_value='Darwin'), \
             patch.object(verifier.subprocess, 'check_output', return_value='Pageouts: 8.\nSwapouts: 4.\n'):
            self.assertEqual(verifier.swap_counters(), {'Swapouts': 4, 'Pageouts': 8})

    def test_machine_has_no_network_identity(self):
        info = executor.machine_info()
        self.assertTrue({'system', 'architecture', 'logical_CPUs', 'python_version'} <= info.keys())
        self.assertFalse({'host', 'IP', 'username', 'password'} & info.keys())


if __name__ == '__main__':
    unittest.main()
