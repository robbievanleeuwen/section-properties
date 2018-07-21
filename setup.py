from setuptools import setup


def description_text():
    return 'A python package for the analysis of arbitrary'\
        ' cross-sections using the finite element method.'


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='sectionproperties',
      version='0.0.1',
      description=description_text(),
      long_description=readme(),
      long_description_content_type='text/markdown',
      classifiers=[
          'Development Status :: 2 - Pre-Alpha',
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
      packages=['sectionproperties'],
      install_requires=[
          'numpy', 'scipy', 'matplotlib'
      ],
      include_package_data=True,
      zip_safe=False)
