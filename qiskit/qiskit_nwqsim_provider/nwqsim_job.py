# -*- coding: utf-8 -*-

# Copyright 2018, IBM.
#
# This source code is licensed under the Apache License, Version 2.0 found in
# the LICENSE.txt file in the root directory of this source tree.

from qiskit.providers import JobV1 as Job
from qiskit.providers import JobError
from qiskit.providers import JobTimeoutError
from qiskit.providers.jobstatus import JobStatus
from qiskit.result import Result

class NWQSimJob(Job):
    _async = False
    def __init__(self, backend, job_id, fn, qobj, options):
        super().__init__(backend, job_id)
        self._result = fn(job_id, qobj, options)
        self._backend = backend
        self._qobj = qobj
        self._options = options
    def submit(self):
        return
    def result(self):
        return self._result
    def status(self):
        return JobStatus.DONE
    def backend(self):
        return self._backend
    def qobj(self):
        return self._qobj
    def options(self):
        return self._options

