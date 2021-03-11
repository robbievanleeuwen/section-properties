import sys
from setuptools import setup
from sectionproperties import __version__ as version


def description_text():
    description = 'A python package for the analysis of arbitrary cross-sections using the'
    description += ' finite element method.'

    return description


def readme():
    with open('README_pypi.rst') as f:
        return f.read()

if sys.version_info < (3, 6):
    sys.exit('Sorry, Python < 3.6 is not supported')

install_requires = ['numpy', 'scipy<1.6', 'matplotlib', 'shapely']

if not (sys.platform == 'win32' or sys.platform == 'cygwin'):
    install_requires.append('pybind11')
    install_requires.append('meshpy')

setup(
    name='sectionproperties',
    version=version,
    description=description_text(),
    long_description=readme(),
    long_description_content_type='text/markdown',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering',
    ],
    url='https://github.com/robbievanleeuwen/section-properties',
    author='Robbie van Leeuwen',
    author_email='robbie.vanleeuwen@gmail.com',
    license='MIT',
    packages=[
        'sectionproperties', 'sectionproperties.analysis', 'sectionproperties.post',
        'sectionproperties.pre', 'sectionproperties.examples', 'sectionproperties.tests'
    ],
    install_requires=install_requires,
    include_package_data=True,
    zip_safe=False
)
