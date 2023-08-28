"""Provider for the NWQSim backend."""

from qiskit.providers import ProviderV1 as Provider
from qiskit.providers.providerutils import filter_backends
from qiskit_nwqsim_provider.nwqsim_simulator import NWQSimSimulator

class NWQSimProvider(Provider):
    """Provider for the NWQSim backend."""
    def __init__(self, token=None):
        super().__init__()
        self.token = token
        self.backends = {'nwqsim': NWQSimSimulator()}

    def backends(self, name=None, **kwargs):
        if name:
            backends = [
                    backend for backend in backends if backend.name() == name]
        return filter_backends(backends, filters=filters, **kwargs)

    def __str__(self):
        return 'NWQSimProvider'
