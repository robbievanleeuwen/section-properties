import sys
from setuptools import setup
from sectionproperties import __version__ as version


def description_text():
    description = 'A python package for the analysis of arbitrary cross-sections using the'
    description += ' finite element method.'

    return description


def readme():
    with open('README.md') as f:
        return f.read()


if sys.version_info[0] < 3 or sys.version_info[0] == 3 and sys.version_info[1] < 7:
    sys.exit('Sorry, Python < 3.7 is not supported')

install_requires = [
    'numpy', 'scipy', 'pytest_check', 'more_itertools', 'matplotlib', 'shapely', 'pybind11',
    'meshpy', 'rhino-shapley-interop~=0.0.1'
]

setup(
    name='sectionproperties',
    version=version,
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
