shewchuk
========

[![](https://dev.azure.com/lycantropos/shewchuk/_apis/build/status/lycantropos.shewchuk?branchName=master)](https://dev.azure.com/lycantropos/shewchuk/_build/latest?branchName=master "Azure Pipelines")
[![](https://codecov.io/gh/lycantropos/shewchuk/branch/master/graph/badge.svg)](https://codecov.io/gh/lycantropos/shewchuk "Codecov")
[![](https://img.shields.io/github/license/lycantropos/shewchuk.svg)](https://github.com/lycantropos/shewchuk/blob/master/LICENSE "License")
[![](https://badge.fury.io/py/shewchuk.svg)](https://badge.fury.io/py/shewchuk "PyPI")

In what follows `python` is an alias for `python3.5` or `pypy3.5`
or any later version (`python3.6`, `pypy3.6` and so on).

Installation
------------

Install the latest `pip` & `setuptools` packages versions
```bash
python -m pip install --upgrade pip setuptools
```

### User

Download and install the latest stable version from `PyPI` repository
```bash
python -m pip install --upgrade shewchuk
```

### Developer

Download the latest version from `GitHub` repository
```bash
git clone https://github.com/lycantropos/shewchuk.git
cd shewchuk
```

Install
```bash
python setup.py install
```

Development
-----------

### Bumping version

#### Preparation

Install
[bump2version](https://github.com/c4urself/bump2version#installation).

#### Pre-release

Choose which version number category to bump following [semver
specification](http://semver.org/).

Test bumping version
```bash
bump2version --dry-run --verbose $CATEGORY
```

where `$CATEGORY` is the target version number category name, possible
values are `patch`/`minor`/`major`.

Bump version
```bash
bump2version --verbose $CATEGORY
```

This will set version to `major.minor.patch-alpha`. 

#### Release

Test bumping version
```bash
bump2version --dry-run --verbose release
```

Bump version
```bash
bump2version --verbose release
```

This will set version to `major.minor.patch`.

### Running tests

Install dependencies
```bash
python -m pip install -r requirements-tests.txt
```

Plain
```bash
pytest
```

Inside `Docker` container:
- with `CPython`
  ```bash
  docker-compose --file docker-compose.cpython.yml up
  ```
- with `PyPy`
  ```bash
  docker-compose --file docker-compose.pypy.yml up
  ```

`Bash` script (e.g. can be used in `Git` hooks):
- with `CPython`
  ```bash
  ./run-tests.sh
  ```
  or
  ```bash
  ./run-tests.sh cpython
  ```

- with `PyPy`
  ```bash
  ./run-tests.sh pypy
  ```

`PowerShell` script (e.g. can be used in `Git` hooks):
- with `CPython`
  ```powershell
  .\run-tests.ps1
  ```
  or
  ```powershell
  .\run-tests.ps1 cpython
  ```
- with `PyPy`
  ```powershell
  .\run-tests.ps1 pypy
  ```
