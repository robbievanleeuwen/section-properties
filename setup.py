import sys
from setuptools import setup


def description_text():
    description = 'A python package for the analysis of arbitrary cross-sections using the'
    description += ' finite element method.'

    return description


def readme():
    with open('README_pypi.rst') as f:
        return f.read()


if sys.version_info[0] < 3 or sys.version_info[0] == 3 and sys.version_info[1] < 6:
    sys.exit('Sorry, Python < 3.6 is not supported')

install_requires = ['numpy', 'scipy<=1.5', 'matplotlib', 'shapely']

if not (sys.platform == 'win32' or sys.platform == 'cygwin'):
    install_requires.append('pybind11')
    install_requires.append('meshpy')

setup(
    name='sectionproperties',
    version='1.0.8',
    description=description_text(),
    long_description=readme(),
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3 :: Only',
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
