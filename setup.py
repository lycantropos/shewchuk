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
        package_data={shewchuk.__name__: ['py.typed']},
        version=shewchuk.__version__,
        description=shewchuk.__doc__,
        long_description=read_file('README.md'),
        long_description_content_type='text/markdown',
        author='Azat Ibrakov',
        author_email='azatibrakov@gmail.com',
        license='MIT License',
        classifiers=[
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: Python :: 3.10',
            'Programming Language :: Python :: 3.11',
            'Programming Language :: Python :: Implementation :: CPython',
            'Programming Language :: Python :: Implementation :: PyPy',
        ],
        url=project_base_url,
        download_url=project_base_url + 'archive/master.zip',
        python_requires='>=3.7'
)
if platform.python_implementation() == 'CPython':
    from glob import glob

    from setuptools import Extension
    from setuptools.command.build_ext import build_ext


    class BuildExt(build_ext):
        def build_extensions(self) -> None:
            compile_args = []
            compiler_type = self.compiler.compiler_type
            if compiler_type == 'unix':
                compile_args.append('-Werror')
            elif compiler_type == 'msvc':
                compile_args.append('/WX')
            for extension in self.extensions:
                extension.extra_compile_args += compile_args
            super().build_extensions()


    parameters.update(
            cmdclass={build_ext.__name__: BuildExt},
            ext_modules=[
                Extension(f'{shewchuk.__name__}._c{shewchuk.__name__}',
                          glob('src/*.c'))
            ],
            zip_safe=False
    )
setup(**parameters)
