# -*- coding: utf-8 -*-

# Copyright 2018, IBM.
#
# This source code is licensed under the Apache License, Version 2.0 found in
# the LICENSE.txt file in the root directory of this source tree.

import functools
import logging
import sys
from concurrent import futures

from qiskit.providers import JobV1 as Job
from qiskit.providers import JobError, JobStatus
#from qiskit.qobj import validate_qobj_against_schema

logger = logging.getLogger(__name__)

def requires_submit(func):
    @functools.wraps(func)
    def _wrapper(self, *args, **kwargs):
        if self._future is None:
            raise JobError("Job is not submitted yet! You have to .submit() first!")
        return func(self, *args, **kwargs)
    return _wrapper

class NWQSimJob(Job):
    if sys.platform in ['darwin', 'win32']:
        _executor = futures.ThreadPoolExecutor()
    else:
        _executor = futures.ProcessPoolExecutor()

    def __init__(self, backend, job_id, fn, qobj):
        super().__init__(backend, job_id)
        self._backend = backend
        self._fn = fn
        self._qobj = qobj
        self._future = None

    def submit(self):
        if self._future is not None:
            raise JobError("The job has already been submitted!")
        #validate_qobj_against_schema(self._qobj)
        self._future = self._executor.submit(self._fn, self._job_id, self._qobj)

    @requires_submit
    def result(self, timeout=None):
        return self._future.result(timeout=timeout)

    @requires_submit
    def get_counts(self, circuit=None, timeout=None):
        return self.result(timeout=timeout).get_counts(circuit)

    @requires_submit
    def cancel(self):
        return self._future.cancel()

    @requires_submit
    def status(self):
        # The order is important here
        if self._future.running():
            _status = JobStatus.RUNNING
        elif self._future.cancelled():
            _status = JobStatus.CANCELLED
        elif self._future.done():
            _status = JobStatus.DONE if self._future.exception() is None else JobStatus.ERROR
        elif self._future._state == 'PENDING':
            _status = JobStatus.QUEUED
        else:
            raise JobError('Unexpected behavior of {0}'.format(
                self.__class__.__name__))
        return _status

    def backend(self):
        return self._backend

    def qobj(self):
        return self._qobj
