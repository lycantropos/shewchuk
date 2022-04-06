import platform
from pathlib import Path

from setuptools import (find_packages,
                        setup)

import shewchuk

project_base_url = 'https://github.com/lycantropos/shewchuk/'


def read_file(path_string: str) -> str:
    return Path(path_string).read_text(encoding='utf-8')


parameters = dict(
        name=shewchuk.__name__,
        packages=find_packages(exclude=('tests', 'tests.*')),
        version=shewchuk.__version__,
        description=shewchuk.__doc__,
        long_description=read_file('README.md'),
        long_description_content_type='text/markdown',
        author='Azat Ibrakov',
        author_email='azatibrakov@gmail.com',
        license='MIT License',
        classifiers=[
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: Python :: 3.10',
            'Programming Language :: Python :: Implementation :: CPython',
            'Programming Language :: Python :: Implementation :: PyPy',
        ],
        url=project_base_url,
        download_url=project_base_url + 'archive/master.zip',
        python_requires='>=3.6')
if platform.python_implementation() == 'CPython':
    from glob import glob
    from typing import Any

    from setuptools import (Command,
                            Extension)
    from setuptools.command.build_ext import build_ext
    from setuptools.command.develop import develop


    class Develop(develop):
        def reinitialize_command(self,
                                 name: str,
                                 reinit_subcommands: int = 0,
                                 **kwargs: Any) -> Command:
            if name == build_ext.__name__:
                kwargs.setdefault('debug', 1)
            result = super().reinitialize_command(name, reinit_subcommands,
                                                  **kwargs)
            if name == build_ext.__name__:
                result.ensure_finalized()
                for extension in result.extensions:
                    extension.undef_macros.append(('NDEBUG',))
            return result


    parameters.update(cmdclass={develop.__name__: Develop},
                      ext_modules=[Extension('_' + shewchuk.__name__,
                                             glob('src/*.c'))],
                      zip_safe=False)
setup(**parameters)
