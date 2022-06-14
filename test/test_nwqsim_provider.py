# -*- coding: utf-8 -*-

# Copyright 2019, IBM.
#
# This source code is licensed under the Apache License, Version 2.0 found in
# the LICENSE.txt file in the root directory of this source tree.

"""Test NWQSim backend."""

from qiskit_nwqsim_provider import NWQSimProvider

from common import QiskitTestCase

class NWQSimProviderTestCase(QiskitTestCase):
    """Tests for the NWQSim Provider."""

    def setUp(self):
        self.provider = NWQSimProvider()
        self.backend_name = 'dmsim'

    def test_backends(self):
        """Test the provider has backends."""
        backends = self.provider.backends()
        self.assertTrue(len(backends) > 0)

    def test_get_backend(self):
        """Test getting a backend from the provider."""
        backend = self.provider.get_backend(name=self.backend_name)
        self.assertEqual(backend.name(), self.backend_name)
