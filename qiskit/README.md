## Install the Python Qiskit interface

1) Setup the Python environment, we need Python>=3.6. It is recommended to use conda for the setup.

2) Install the importlib_resources package for locating the nwqsim lib. 
```bash
pip install importlib_resources
```

3) Build the package
```bash
python setup build
```

4) Install the package
```bash
python setup install
```

5) Test the setup
```bash
python ghz.py
```

Please find how to call the simulator in the examplar ghz.py.
